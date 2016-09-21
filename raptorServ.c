/* -*- c-file-style: "Ellemtel" -*- */
/* Copyright (C) 2015   Canada-France-Hawaii Telescope Corp.          */
/* This program is distributed WITHOUT any warranty, and is under the */
/* terms of the GNU General Public License, see the file COPYING      */
/*!**********************************************************************
 * $Id$
 * $Locker$
 *
 * DESCRIPTION
 *
 *    Server process which connects to the Raptor Camera via a 
 *    CameraLink interface and enables FITS files to be sent to STDOUT
 *    These FITS files can be intercepted and transferred to a FITS 
 *    server where they can be used by downstream processes.
 *
 *    The camera get/set functions were authored by Chi-Hung Yang.
 *    They were extracted from an sgc.c program and modified to 
 *    incorporate CFHT logging standards.
 *
 *    Parts of the code related to ISU (Image Stabilization Unit) control
 *    were authored by Marie Larrieu. They call functions from the libisu
 *    library, which is a PowerDAQ interface to the ISU steering mechanism.
 * HISTORY
 *
 * $Log$
 *
 *********************************************************************!*/
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <netdb.h>
#include <sys/socket.h>

#include "cli/cli.h"
#include "fh/fh.h"
#include "fh/fh_registry.h"
#include "sockio/sockserv.h"
#include "ss/linked_list.h"
#include "ssapi/ss_api.h"
#include "ssapi/ss_error.h"

#include "edtinc.h"
#include "pciload.h"
#include "clsim_lib.h"

#include "mpfit/mpfit.h"

#include "isu.h"

#define RAPTOR_PORT "915"  /* Port name or number on which to listen */
#define RAPTOR_CONFIG "/cfht/conf/raptor.conf"


/*
 * Telescope pointing information
 */
#define SS_TELESCOPE_RA "/t/status/currentRA"
#define SS_TELESCOPE_DEC "/t/status/currentDEC"
#define SS_TELESCOPE_EQ "/t/status/currentEQ"

/*
 * Serial channel interface parameters
 */
#define SERIALTIMEOUT 6
#define SERBUFSIZE 512
#define UNIT       0
#define BAUD       115200
#define CHANNEL    0
#define OFF        0
#define ON         1
#define LOWGAIN    0
#define HIGHGAIN   1

#define USER_TIMEOUT 20000

#define DEFAULT_TEC_SETPOINT -40.0
#define DEFAULT_DIGITAL_GAIN 1
#define DEFAULT_EXPOSURE_TIME 10.0
#define DEFAULT_FRAME_RATE 50.0
#define MAX_SAVE_COUNT 1000000
#define SATURATION 150000

#define QUIT_CMD "QUIT"
#define BYE_CMD "BYE"
#define EXIT_CMD "EXIT"
#define LOGOUT_CMD "LOGOUT"
#define SHUTDOWN_CMD "SHUTDOWN"
#define FRAMERATE_CMD "FRAMERATE"
#define EXPTIME_CMD "EXPTIME"
#define TEC_CMD "TEC"
#define TEMP_CMD "TEMP"
#define ROI_CMD "ROI"
#define VIDEO_CMD "VIDEO"
#define SAVE_CMD "SAVE"
#define GUIDE_CMD "GUIDE"
#define ISU_CMD "ISU"
#define PASS_CHAR '.'
#define FAIL_CHAR '!'
#define OOB_CHAR '*'
#define PERR_CHAR '?'


#define GUIDE_SIZE_X 32
#define GUIDE_SIZE_Y 32

#define PIXSCALE 0.128
//int   fwhm_flag=0;
//float fwhm_x, fwhm_y;


#define SOCKSERV_POLL_INTERVAL 1 /* Corresponds to 1/100 second */

/*
 * BEWARE THAT
 * if DEBUG is defined, some printf are done on a file and on stderr
 * so this is memory and time consuming
 */
//#define DEBUG

/*
 * BEWARE THAT
 * if SIM_STAR is defined, the centroid algorithm is not called 
 * and a software gaussian star simulation is sent to the ISU 
 * This definition should be COMMENTED FOR NORMAL OPERATIONS !
 */
#define SIM_STAR

/*
 * Per-client information.  Multiple clients can stay connected to the
 * server at once, and the server must keep track of the current exposure
 * type, exposure time, and raster settings for each of them since any
 * of them may be next to send a "go" (and then the last-requested settings
 * on *that* socket need to be re-sent to libfli.)
 */
typedef struct
{
      char *hostname;
      unsigned char remote_ip[4];
      time_t connect_ts;
} client_info_t;


/*
 * Structure used to specify server specific information.
 */
typedef struct {
      linked_list *client_list;
      sockserv_t *raptor_serv;
      int serv_done;
      Dependent *dd_p;
      EdtDev *edt_p;
      EdtDev *pdv_p;
      float frame_rate;  /* Hz */
      float exposure_time;
      float tec_setpoint;
      float temp;
      int image_width;
      int image_height;
      int win_x0;
      int win_y0;
      float guide_xoff;
      float guide_yoff;
      BOOLEAN video_on;
      BOOLEAN isu_on;
      BOOLEAN guide_on;
      char *fits_comment;
      int fwhm_flag;
      float fwhm_x;
      float fwhm_y;
      int frame_sequence;
      int frame_save_count;
} server_info_t;


/*
 * Server information structure instance
 */
static server_info_t *serv_info;

/*MPFIT STRUCTURE - DEFINE THE PRIVATE STRUCTURE FOR THE DATA, ERRORS,
  COORDINATES ETC. ANY 2D BEHAVIOUR IS HERE.  REMOVED X AND Y BECAUSE WE
  CAN GET THESE, ASSUMING IMAGE IS SQUARE */

struct vars_struct {
  double *flux;     //DATA
  double *ferr;     //ESTIMATE OF ERROR
};
//--------------------------------------------------//

/*
 * This function calls set_analog_slope in a threaded fashion
 *
 */
static void *set_analog_slope_thread(void* p_args)
{
   if (p_args != NULL)
   {
      thdata *data;
      data = (thdata *) p_args;
      if (set_analog_slope(data->arg1,
                           data->arg2, data->arg3,
                           data->arg4, data->arg5) == FAIL) {
         cfht_logv(CFHT_MAIN, CFHT_WARN,
                   "(%s:%d) %s: Failed setting analog output slope with "
                   "arguments  %lf %lf %lf %lf %lf",
                 __FILE__, __LINE__, __FUNCTION__, data->arg1,
                  data->arg2, data->arg3, data->arg4, data->arg5);
      }
   }
   else {
     cfht_logv(CFHT_MAIN, CFHT_ERROR,
               "(%s:%d) %s: Missing arguments",
                 __FILE__, __LINE__, __FUNCTION__);

   }
   pthread_exit(NULL);
}

/* QUICK ROUTINE TO PRINT THE RESULTS FOR A PARTICULAR MPFIT - LABEL, NUMBER OF ITERATIONS
   AND THE FINAL VALUE OF THE PARAMETERS. FOR DEBUGGING ONLY */

void printresult(double *x, mp_result *result)
{
  int i;

  if ((x == 0) || (result == 0)) return;  //NO RESULTS

  fprintf(stderr,"%d ",result->niter);
  for (i=0; i<result->npar; i++)
    {
	  fprintf(stderr," %lf ",x[i]);
    }
  fprintf(stderr,"\n");


}

#ifdef SIM_STAR
/*
 * This polar form of the Box-Muller transformation
 * We start with two independent random numbers, x1 and x2, 
 * which come from a uniform distribution (in the range from 0 to 1).
 * Then apply the above transformations to get two new independent 
 * random numbers which have a Gaussian distribution with zero mean 
 * and a standard deviation of one. 
 * Source: http://www.design.caltech.edu/erik/Misc/Gaussian.html
 */
void sampleNormal(float* y1, float* y2)
{
   float x1, x2, w;
 
   do {
      x1 = 2.0 * ((double)rand() / (RAND_MAX)) - 1.0;
      x2 = 2.0 * ((double)rand() / (RAND_MAX)) - 1.0;
      w = x1 * x1 + x2 * x2;
   } while ( w >= 1.0 );

   w = sqrt( (-2.0 * log( w ) ) / w );
   *y1 = x1 * w;
   *y2 = x2 * w;
}
#endif

//--------------------------------------------------//
static PASSFAIL
gaussfunc2d(int m, int n, double *p, double *dy, double **dvec, void *vars)
{

  int i,j;    //COUNTERS, SIZE OF IMAGE
  int n1;     //SIZE OF IMAGE, ASSUMING A SQUARE
  struct vars_struct *v = (struct vars_struct *) vars;  //STRUCTURE TO PASS PRIVATE DATA
  double *flux;  //DATA VALUES
  double *ferr;  //ERRORS
  double xc, yc;        //THE VALUE OF THE CENTRE POINTS OF THE IMAGE

  //SET THE LOCAL VARIABLES TO THE STRUCTURE.

  flux=v->flux;
  ferr = v->ferr;

  //ASSUMING A SQUARE IMAGE, GET THE DIMENSIONS OF EACH SIDE

  n1=(int)sqrt(m);


  //CYCLE THROUGH THE VALUES. THE DATA/RESIDUALS ARE ONE D,
  //MAP THE COORDINATES TO THAT ASSUMING A SQUARE

  for (i=0;i<n1;i++)
    {
      for (j=0;j<n1;j++)
    {
	  //CENTER VALUES
	  xc=i-p[0];
	  yc=j-p[1];

	  //EQUATION ASSUMING INDEPENDENT FWHM IN X AND Y DIRECTIONS

	  dy[i*n1+j] = ((flux[i*n1+j] - p[4]*exp(-0.5*(xc*xc/(p[2]*p[2]*0.180337)+yc*yc/(p[3]*p[3]*0.180337)))-p[5]));
    }
    }
  return PASS;
}


#define ELEM_SWAP(a,b) {register float t=(a);(a)=(b);(b)=t; }

/** \fn  float GetF(float arr[], int n)
 * \brief This function is the quick_select routine based on the algorithm found
 * in Numerical Recipes in C.
 *
 * \param arr[] is the image
 * \param n is the number of pixels in the image
 * \return median value of an array
 */
double GetMedian(double arr[], int n){
	int low, high;
	int median;
	int middle, ll, hh;

	low = 0; high = n-1; median = (low + high) / 2;
	for (;;) {

		if (high <= low) /* One element only */
			return arr[median] ;

		if (high == low + 1) { /* Two elements only */
			if (arr[low] > arr[high])
				ELEM_SWAP(arr[low], arr[high]) ;
			return arr[median] ;
		}

		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high])
			ELEM_SWAP(arr[middle], arr[high])
		if (arr[low] > arr[high])
			ELEM_SWAP(arr[low], arr[high])
		if (arr[middle] > arr[low])
			ELEM_SWAP(arr[middle], arr[low])

		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP(arr[middle], arr[low+1]) ;

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do ll++; while (arr[low] > arr[ll]) ;
			do hh--; while (arr[hh] > arr[low]) ;

				if (hh < ll) break;

				ELEM_SWAP(arr[ll], arr[hh])
		}

		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP(arr[low], arr[hh])

		/* Re-set active partition */
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1;
	}

}
#undef ELEM_SWAP
/*
 * Simple centroid calculation on the image
 */
static void
calculateCentroid(unsigned short *image, int columns, int rows,
		  float *xc, float *yc) {

   int i;
   int j;
   int val;
   float sum,median;
   double arr[columns*rows];

   //arr=(double *)malloc(columns*rows*sizeof(double));

   for (i=0;i<columns*rows;i++){
	   arr[i]=image[i];
   }
   median=GetMedian(arr,columns*rows);
   //fprintf(stderr," median=%f ",median);

   *xc = 0;
   *yc = 0;
   sum = 0;
   for (i = 0; i < rows; i++) {
      for (j = 0; j < columns; j++) {
	 val = (image)[i * columns + j]-median;
	 if (val > 0) {
	    *xc += j * val;
	    *yc += i * val;
	    sum += val;
	 }
      }
   }
   if (sum > 0) {
      *xc /= sum;
      *yc /= sum;
   }
   else {
      *xc = columns / 2.0;
      *yc = rows / 2.0;
   }

   //free(arr);
}


/*
 * MPFIS method for centroid calculation on the image
 */
void *calculateCentroidMPFIT(unsigned short *image, int columns, int rows,
		  float *xc, float *yc) {

	struct vars_struct v; //PRIVATE STRUCTURE WITH DATA/FUNCTION INFORMATION
	mp_result result;     //STRUCTURE WITH RESULTS
	mp_par pars[6];	//VARIABLE THAT HOLDS INFORMATION ABOUT FIXING PARAMETERS - DATA TYPE IN MPFIT

	double median;

	double *ferr=NULL;

	// Array for median calculation
	double arr[columns*rows];

	double perror[6];	//ERRORS IN RETURNED PARAMETERS

	int status;          //EXIT STATUS
	int i,j,k=0;
	int npoints,hpoint;

	/*
	 *   First step, estimate the center of the point using Center of Mass
	 */
	float xest = 0;
	float yest = 0;
	calculateCentroid(image, GUIDE_SIZE_X, GUIDE_SIZE_Y,&xest, &yest);
	//fprintf(stderr,"x=%f y=%f ",xest,yest);

	/*
	 *  Cut out the region near the point
	 */
	int fpix[2];                       //FIRST PIXELS OF SUBREGION [X,Y]
	int lpix[2];                       //LAST PIXELS OF SUBREGION [X,Y]
	int subx, suby,np;
	double *subimage;

	fpix[0]=xest-GUIDE_SIZE_X/4;
	fpix[1]=yest-GUIDE_SIZE_Y/4;
	lpix[0]=xest+GUIDE_SIZE_X/4-1;
	lpix[1]=yest+GUIDE_SIZE_Y/4-1;

	if (xest-GUIDE_SIZE_X/4 < 0) fpix[0]=0;
	if (yest-GUIDE_SIZE_X/4 < 0) fpix[1]=0;
	if (xest-GUIDE_SIZE_X/4 > 32) lpix[0]=32;
	if (yest-GUIDE_SIZE_X/4 > 32) lpix[1]=32;


	//GET THE DIMENSIONS OF THE SUBREGION
	subx=lpix[0]-fpix[0]+1;
	suby=lpix[1]-fpix[1]+1;

	// Copy the central region
	subimage = malloc(subx*suby*sizeof(double));
	ferr = malloc(subx*suby*sizeof(double));
	//fprintf(stderr,"subx=%i %i \n",subx,suby);

	for (i=fpix[0];i<fpix[0]+subx;i++){
		for (j=fpix[1];j<fpix[1]+suby;j++){
			subimage[k]=(double)image[j*columns+i];
			ferr[k]=1.0;
			k++;
		}
	}

	npoints=columns*rows;
	np=subx*suby;
	hpoint=columns/2;

	//ferr = malloc(npoints*sizeof(double));
	//arr = malloc(npoints*sizeof(double));
	for (i=0;i<npoints;i++){
		//ferr[i]=1.0;
		arr[i]=(double)image[i];
	}
	median=GetMedian(arr,columns*rows);

	double p[] = {xest-fpix[0],yest-fpix[1],2.5,2.5,12800.0,median};

	memset(&result,0,sizeof(result));
	result.xerror = perror;
	memset(pars,0,sizeof(pars));
	//fprintf(stderr,"init=%f %f\n",xest-fpix[0],yest-fpix[1]);

	v.ferr = ferr;
	v.flux = subimage;

   	//pars[1].fixed = 0;
	pars[2].fixed = 1;
	pars[3].fixed = 1;
	//pars[4].fixed = 1;
	//pars[5].fixed = 1;
	pars[5].fixed = 1;



	status= mpfit(gaussfunc2d, np, 6, p, pars, 0, (void *) &v, &result);


	//fprintf(stderr,"total time=%f \n",last_ts-current_ts);

	//printresult(p, &result);

	free(subimage);
	free(ferr);

	if (fpix[1]+p[1] < 0){
		*yc=yest;
	} else {
		*yc=fpix[1]+p[1];
	}

	if (fpix[0]+p[0] < 0){
		*xc=xest;
	} else {
		*xc=fpix[0]+p[0];
       	}

}

/*
 * This function is used to calculate the FWHM of the stellar point.
 */
void *calculatePointFWHM(unsigned short *image, int columns, int rows)
{

	struct vars_struct v; //PRIVATE STRUCTURE WITH DATA/FUNCTION INFORMATION
	mp_result result;     //STRUCTURE WITH RESULTS
	mp_par pars[6];	//VARIABLE THAT HOLDS INFORMATION ABOUT FIXING PARAMETERS - DATA TYPE IN MPFIT

	double median;

	double *ferr=NULL;

	// Array for median calculation
	double arr[columns*rows];

	double perror[6];	//ERRORS IN RETURNED PARAMETERS

	int status;          //EXIT STATUS
	int i,j,k=0;
	int npoints,hpoint;

	/*
	 *   First step, estimate the center of the point using Center of Mass
	 */
	float xest = 0;
	float yest = 0;
	calculateCentroid(image, GUIDE_SIZE_X, GUIDE_SIZE_Y,&xest, &yest);

	#ifdef DEBUG
	fprintf(stderr,"x=%f y=%f \n",xest,yest);
	#endif
	/*
	 *  Cut out the region near the point
	 */
	int fpix[2];                       //FIRST PIXELS OF SUBREGION [X,Y]
	int lpix[2];                       //LAST PIXELS OF SUBREGION [X,Y]
	int subx, suby,np;
	double *subimage;

	fpix[0]=xest-GUIDE_SIZE_X/4;
	fpix[1]=yest-GUIDE_SIZE_Y/4;
	lpix[0]=xest+GUIDE_SIZE_X/4-1;
	lpix[1]=yest+GUIDE_SIZE_Y/4-1;

	if (xest-GUIDE_SIZE_X/4 < 0) fpix[0]=0;
	if (yest-GUIDE_SIZE_X/4 < 0) fpix[1]=0;
	if (xest-GUIDE_SIZE_X/4 > 32) lpix[0]=32;
	if (yest-GUIDE_SIZE_X/4 > 32) lpix[1]=32;


	//GET THE DIMENSIONS OF THE SUBREGION
	subx=lpix[0]-fpix[0]+1;
	suby=lpix[1]-fpix[1]+1;

	// Copy the central region
	subimage = malloc(subx*suby*sizeof(double));
	ferr = malloc(subx*suby*sizeof(double));
	//fprintf(stderr,"subx=%i %i \n",subx,suby);

	for (i=fpix[0];i<fpix[0]+subx;i++){
		for (j=fpix[1];j<fpix[1]+suby;j++){
			subimage[k]=(double)image[j*columns+i];
			ferr[k]=1.0;
			k++;
		}
	}

	npoints=columns*rows;
	np=subx*suby;
	hpoint=columns/2;

	//ferr = malloc(npoints*sizeof(double));
	//arr = malloc(npoints*sizeof(double));
	for (i=0;i<npoints;i++){
		//ferr[i]=1.0;
		arr[i]=(double)image[i];
	}
	median=GetMedian(arr,columns*rows);

	double p[] = {xest-fpix[0],yest-fpix[1],2.5,2.5,12800.0,median};

	memset(&result,0,sizeof(result));
	result.xerror = perror;
	memset(pars,0,sizeof(pars));
	//fprintf(stderr,"init=%f %f\n",xest-fpix[0],yest-fpix[1]);

	v.ferr = ferr;
	v.flux = subimage;

   	//pars[1].fixed = 0;
	//pars[2].fixed = 1;
	//pars[3].fixed = 1;
	//pars[4].fixed = 1;
	//pars[5].fixed = 1;
	pars[5].fixed = 1;



	status= mpfit(gaussfunc2d, np, 6, p, pars, 0, (void *) &v, &result);


	//fprintf(stderr,"total time=%f \n",last_ts-current_ts);

	//printresult(p, &result);
	serv_info->fwhm_x=p[2];
	serv_info->fwhm_y=p[3];

	free(subimage);
	free(ferr);

}



/* 
 * Advance past leading whitespace in a string 
 */
char *
ltrim(char *str) 
{
   
   char *p;
   
   if (str == NULL) 
      return NULL;
   
   for (p = str; (*p && isspace(*p)); p++)
      ;

   return p;
}


/* 
 * Trim off trailing whitespace in a string 
 */
char *
rtrim(char *str) 
{

   int i;

   if (str != NULL) {
      for (i = strlen(str) - 1; (i >= 0 && isspace(str[i])); i--)
	 ;
      str[++i] = '\0';
   }   

   return str;
}


/* 
 * Trim off all leading and trailing whitespace from a string 
 */
char *
trim(char *str) 
{

   return rtrim(ltrim(str));
}


char ** 
stringSplit(char *a_str, const char a_delim)
{
   char** result    = 0;
   size_t count     = 0;
   char* tmp        = a_str;
   char* last_comma = 0;
   char delim[2];
   delim[0] = a_delim;
   delim[1] = 0;

   //printf("%s\n",a_str);
   /* Count how many elements will be extracted. */
   while (*tmp)
   {
      if (a_delim == *tmp)
      {
	 count++;
	 last_comma = tmp;
      }
      tmp++;
   }
   
   /* Add space for trailing token. */
   count += last_comma < (a_str + strlen(a_str) - 1);
   
   /* Add space for terminating null string so caller
      knows where the list of returned strings ends. */
   count++;
   
   result = malloc(sizeof(char*) * count);
   
   if (result)
   {
      size_t idx  = 0;
      char* token = strtok(a_str, delim);
      
      while (token)
      {
	 assert(idx < count);
	 *(result + idx++) = strdup(token);
	 token = strtok(0, delim);
      }
      assert(idx == count - 1);
      *(result + idx) = 0;
   }
   
   return result;
}


/*
 * Deallocate a client object
 */
static void 
deallocateClient(void *data) 
{
   client_info_t *cinfo = (client_info_t *)data;

   if (data == NULL) {
      return;
   }

   if (cinfo->hostname != NULL) {
      free(cinfo->hostname);
   }
   free(cinfo);
}


/*
 * Callback function used to compare two pointers.  This function is 
 * used within the list comparison utilities.
 */
static int
comparePointer(const void *a, const void *b)
{
   if (a < b)
      return -1;
   else if (a > b)
      return 1;
   else
      return 0;
}


/*
 * Perform a write across the serial channel and read back the response
 */
static PASSFAIL 
pdvSerialWriteRead(char *ibuf_p, char **response){

   int timeout = SERIALTIMEOUT;
   int ret, i;
   int length=0;
   static char temp[256];
   char buf[SERBUFSIZE+1];
   u_char hbuf[SERBUFSIZE];
   u_char lastbyte, waitc;
   u_int val;
  

   temp[0] = '\0';

   /* 
    * open a handle to the device     
    */
   if (serv_info->pdv_p == NULL) {
      if ((serv_info->pdv_p = pdv_open_channel(EDT_INTERFACE, 
					       UNIT, CHANNEL)) == NULL) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) pdv_open_channel request failed",
		   __FILE__, __LINE__);
	 return FAIL;
      }
   }
   
   /*
    * Enable reading across the serial interface
    */
   pdv_serial_read_enable(serv_info->pdv_p);
   
   /*
    * Get the timeout value from EDT card configuration if it is too short
    */
   if (timeout < 1) {
      timeout = serv_info->pdv_p->dd_p->serial_timeout;
   }
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) serial timeout value = %d",
	     __FILE__, __LINE__, timeout);

   /*
    * Set the baud rate for the serial channel
    */
   pdv_set_baud(serv_info->pdv_p, BAUD);
   
   /* 
    * Flush any junk on the interface 
    */
   (void)pdv_serial_read(serv_info->pdv_p, buf, SERBUFSIZE);
   
   i = 0;
   
   strip_newline(ibuf_p);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) serial command request = %s",
	     __FILE__, __LINE__, ibuf_p);

   /*
    * Process the serial input request in order to detect possible errors
    */ 
   while (*ibuf_p) {
      while ((*ibuf_p == ' ') || (*ibuf_p == '\t')) {
	 ++ibuf_p;
      }
      
      if (*ibuf_p == '\0') {
	 break;
      }
      
      if (sscanf(ibuf_p, "%x", &val) != 1) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) error reading input byte %d", 
		   __FILE__, __LINE__, i);
	 i = 0;
	 break;
      }
      
      if (val > 0xff) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) hex string format error -- expect hex bytes"
		   " separated by spaces, e.g. '00 a0 ff ...",
		   __FILE__, __LINE__);
	 i = 0;
	 break;
      }
      else {
	 hbuf[i++] = val;
      }
      
      while ((*ibuf_p != ' ') && (*ibuf_p != '\t') && (*ibuf_p != '\0')) {
	 ++ibuf_p;
      }
   }


   /*
    * using pdv_serial_binary_command instead of
    * pdv_serial_write because it prepends a 'c' if FOI
    */
   if (i) {
      if (pdv_serial_binary_command(serv_info->pdv_p, (char *)hbuf, i) != 0) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) can not send serial binary "
		   "command to camera", __FILE__, __LINE__);
	 return FAIL;
      }
   }
   
   /*
    * serial_timeout comes from the config file (or -t override flag in
    * this app), or if not present defaults to 500 unless readonly
    * defaults to 60000
    */
   pdv_serial_wait(serv_info->pdv_p, timeout, 64);
   
   /*
    * Handle the response
    */
   lastbyte = 0;
   do {
      ret = pdv_serial_read(serv_info->pdv_p, buf, SERBUFSIZE);
      cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		"(%s:%d) read returned %d",
		__FILE__, __LINE__, ret);

      if (*buf) {
	 lastbyte = (u_char)buf[strlen(buf)-1];
      }
      
      if (ret != 0) {
	 int i;

	 buf[ret + 1] = 0;
	 
	 if (ret) {
	    for (i = 0; i < ret; i++) {
	       sprintf(temp, "%s %02x", temp,(u_char) buf[i]);
	    }
	 }
	 length += ret;
      }
      
      if (serv_info->pdv_p->devid == PDVFOI_ID) {
	 ret = pdv_serial_wait(serv_info->pdv_p, 500, 0);
      } else if (pdv_get_waitchar(serv_info->pdv_p, &waitc) && 
		 (lastbyte == waitc)) {
	 ret = 0; /* jump out if waitchar is enabled/received */
      }
      else {
	 ret = pdv_serial_wait(serv_info->pdv_p, 500, 64);
      }
      
   } while (ret > 0);
   
//   pdv_close(serv_info->pdv_p);
//   serv_info->pdv_p = NULL;
   
   *response = trim(temp);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) response = %s", __FILE__, __LINE__, *response);
   
   return PASS;
}


/*
 * Check the camera status
 */
static PASSFAIL 
checkCameraStatus(char **response){

   char send_string[256];

   /*
    * Send the "Get system status" command to the camera
    */
   snprintf(send_string, sizeof(send_string)-1, "49 50 19");
   if (pdvSerialWriteRead(send_string, response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) there is no response from the camera...check power",
		__FILE__, __LINE__);
      return FAIL;
      
   }
   
   /*
    * Set the "set system status (=0x53)" command to the camera
    */
   snprintf(send_string, sizeof(send_string)-1, "4f 53 50 4c");
   if (pdvSerialWriteRead(send_string, response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) there is no response from the camera...check power",
		__FILE__, __LINE__);
      return FAIL;
      
   }

   /*
    * Check the response against the expected response
    */
   if (strcmp ("50 4c", *response) != 0) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) expected response not received from the camera",
		__FILE__, __LINE__);
      return FAIL;
   }

   return PASS;
}


/*
 * Set the non-uniform image correction state
 */
static PASSFAIL
setGuiderNUC(int mode){

   char send_string[256];
   char *response;

   if (mode == 0){
      snprintf(send_string, sizeof(send_string)-1, "53 e0 02 f9 01 50 19");
      if (pdvSerialWriteRead(send_string, &response) != PASS) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) no response from camera when setting NUC state", 
		   __FILE__, __LINE__);
	 return FAIL;
      }
      
      /* 
       * Check the response against the expected response
       */
      if (strcmp ("50 19", response) != 0){
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		   "(%s:%d) not able to turn off camera non-uniform"
		   " correction", __FILE__, __LINE__);
	 return FAIL;
      }
   } else {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) not able to turn on camera non-uniform correction."
		"  Functionality not implemented yet", __FILE__, __LINE__);
      return FAIL;
   }
   
   return PASS;
}


/*
 * Turn on or off the guider auto level
 */
static PASSFAIL
setGuiderAutoLevel(int mode) {

   char send_string[256];
   char *response;
   
   if (mode == 0) {
      snprintf(send_string, sizeof(send_string)-1, "53 e0 02 23 00 50 c2");
      if (pdvSerialWriteRead(send_string, &response) != PASS) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) no response from camera when setting auto level"
		   " state", __FILE__, __LINE__);
	 return FAIL;
      }
      
      if (strcmp ("50 c2", response) != 0) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		   "(%s:%d) not able to turn off camera automatic level", 
		   __FILE__, __LINE__);
	 return FAIL;
      }
   } else {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) not able to turn on camera automatic level."
		"  Functionality not implemented yet", __FILE__, __LINE__);
      return FAIL;
   }

   return PASS;
}


/*
 * Enable the TEC cooler
 */
static PASSFAIL
enableGuiderTEC(void) {

   char send_string[256];
   char *response;
  
   snprintf(send_string, sizeof(send_string)-1, "53 e0 02 00 81 50 60");
   if (pdvSerialWriteRead(send_string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) no response from camera when enabling the TEC cooler",
		__FILE__, __LINE__);
      return FAIL;
   }

   if (strcmp ("50 60", response) != 0) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) can not enable TEC cooler", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   return PASS;
}


/*
 * Set the guider TEC temperature set point
 */
static PASSFAIL
setGuiderTECPoint(float temp){

   char hexstring[50];
   char longreps[256];
   char strtmp[3];
   
   char *response;
   char string[256];
   char cmdstring[256];
   char **tokens;
   
   unsigned int dac40d, dac0d;
   unsigned int comd0[6]={0x53,0xe0,0x02,0xfb,0x00,0x50};
   unsigned int comd1[6]={0x53,0xe0,0x02,0xfa,0x00,0x50};
   
   unsigned int valtmp, csum0, csum1;
   
   
   float co,slope,count;
   
   hexstring[0] = '\0';

   /* 
    * Getting manufacturing data from camera 
    */
   snprintf(string, sizeof(string)-1, "53 ae 05 01 00 00 02 00 50 ab");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 af 12 50 be");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s", 
	    *(tokens+17), *(tokens+16));
   dac40d = strtoul(hexstring, 0, 16);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) DAC40d --> HEX = %s value=%i",
	     __FILE__, __LINE__, hexstring, dac40d);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s",
	    *(tokens+15), *(tokens+14));
   dac0d = strtoul(hexstring, 0, 16);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) DAC0d --> HEX = %s value=%i",
	     __FILE__, __LINE__, hexstring, dac0d);
   free(tokens);
   
   slope = (dac40d - dac0d) / (40);
   co = dac0d;
   count = temp * slope + co;
   
   snprintf(hexstring, sizeof(hexstring)-1, "%04lx", (unsigned long)(count));
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) hexstring = %s value=%f",
	     __FILE__, __LINE__, hexstring, count);
   
   /*
    * Preparing serial command 
    */

   /* 
    * Command line 1 
    */
   strncpy(strtmp, hexstring, 2);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);
   csum0 = comd0[0] ^ comd0[1] ^ comd0[2] ^ comd0[3] ^ valtmp ^ comd0[5];
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) comd= %02x %02x %02x %02x %02x %02x %02x",
	     __FILE__, __LINE__, comd0[0], comd0[1], comd0[2], comd0[3],
	     valtmp, comd0[5], csum0);
   
   snprintf(cmdstring, sizeof(cmdstring)-1, 
	    "%02x %02x %02x %02x %02x %02x %02x",
	    comd0[0], comd0[1], comd0[2], comd0[3], valtmp, comd0[5], csum0);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) command=%s", __FILE__, __LINE__, cmdstring);
   
   if (pdvSerialWriteRead(cmdstring, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   /* 
    * Command line 2 
    */
   snprintf(hexstring, sizeof(hexstring), "%04lx", (unsigned long)(count));
   strncpy(strtmp, hexstring + 2, 3);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);
   
   csum1 = comd1[0] ^ comd1[1] ^ comd1[2] ^ comd1[3] ^ valtmp ^ comd1[5];
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) comd= %02x %02x %02x %02x %02x %02x %02x",
	     __FILE__, __LINE__, comd1[0], comd1[1], comd1[2], comd1[3],
	     valtmp, comd1[5], csum1);
   
   snprintf(cmdstring, sizeof(cmdstring)-1, 
	    "%02x %02x %02x %02x %02x %02x %02x",
	    comd1[0], comd1[1], comd1[2], comd1[3], valtmp, comd1[5], csum1);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) command=%s", __FILE__, __LINE__, cmdstring);
   
   if (pdvSerialWriteRead(cmdstring, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
      
   return PASS;
}


/* 
 * Get the guider TEC set point
 */
static PASSFAIL 
getGuiderTECPoint(float *temp) {

   char hexstring[50];
   char longreps[256];
   
   char *response;
   char string[256];
   char **tokens;
   unsigned long dac40d, dac0d;
   
   float value, slope, co;

   hexstring[0] = '\0';
   
   /* 
    * Getting manufacturing data from camera 
    */
   snprintf(string, sizeof(string)-1, "53 ae 05 01 00 00 02 00 50 ab");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 af 12 50 be");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s", 
	    *(tokens+17), *(tokens+16));
   dac40d = strtoul(hexstring, 0, 16);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) DAC40d --> HEX = %s value=%lu",
	     __FILE__, __LINE__, hexstring, dac40d);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s",
	    *(tokens+15), *(tokens+14));
   dac0d = strtoul(hexstring, 0, 16);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) DAC0d --> HEX = %s value=%lu",
	     __FILE__, __LINE__, hexstring, dac0d);
   free(tokens);
   
   /* Reading Current TEC setpoint */
   snprintf(string, sizeof(string)-1, "53 e0 01 fb 50 19");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(hexstring, sizeof(hexstring)-1, "%s", *(tokens));
   free(tokens);
   
   snprintf(string, sizeof(string)-1, "53 e0 01 fa 50 18");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(string, sizeof(string)-1, "%s", hexstring);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s", string, *(tokens));
   free(tokens);
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) HEX = %s", __FILE__, __LINE__, hexstring);
   value = (float)strtoul(hexstring, 0, 16);
   slope = (40.0) / ((float)dac40d - (float)dac0d);
   co = -(slope * (float)dac0d);
   
   *temp = (slope * value) + co;
   
   return PASS;
}


/*
 * Set the frame rate for the guider
 */
static PASSFAIL
setGuiderFrameRate(double count) {

   unsigned char hexstring[20];
   char string[80];
   char *response;
   
   char strtmp[3];
   
   unsigned int comd0[6]={0x53,0xe0,0x02,0xdd,0x00,0x50};
   unsigned int comd1[6]={0x53,0xe0,0x02,0xde,0x00,0x50};
   unsigned int comd2[6]={0x53,0xe0,0x02,0xdf,0x00,0x50};
   unsigned int comd3[6]={0x53,0xe0,0x02,0xe0,0x00,0x50};
   
   unsigned int valtmp, csum0, csum1, csum2, csum3;
   
   snprintf(hexstring, sizeof(hexstring)-1, "%08lx",
	    (unsigned long)(40e8/(unsigned long)(count*100)));
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) hexstring = %s value=%f count = %lu", 
	     __FILE__, __LINE__, hexstring, count,
	     (unsigned long)(40e8/(unsigned long)(count*100)));
   
   /*
    * Preparing serial command 
    */

   /* 
    * Command line 1 
    */
   strncpy(strtmp, hexstring, 2);
   strtmp[2] = '\0';

   valtmp = strtoul(strtmp, 0, 16);

   csum0 = comd0[0] ^ comd0[1] ^ comd0[2] ^ comd0[3] ^ valtmp ^ comd0[5];

   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) comd= %02x %02x %02x %02x %02x %02x %02x",
	     __FILE__, __LINE__, comd0[0], comd0[1], comd0[2],
	     comd0[3], valtmp, comd0[5], csum0);
   
   snprintf(string, sizeof(string)-1, "%02x %02x %02x %02x %02x %02x %02x",
	    comd0[0], comd0[1], comd0[2], comd0[3], valtmp, comd0[5], csum0);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) command = %s", __FILE__, __LINE__, string);
   
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   /* 
    * Command line 2 
    */
   snprintf(hexstring, sizeof(hexstring)-1, "%08lx",
	    (unsigned long)(40e8/(unsigned long)(count*100)));
   strncpy(strtmp, hexstring + 2, 3);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);
   
   csum1 = comd1[0] ^ comd1[1] ^ comd1[2] ^ comd1[3] ^ valtmp ^ comd1[5];

   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) comd= %02x %02x %02x %02x %02x %02x %02x",
	     __FILE__, __LINE__, comd1[0], comd1[1], comd1[2],
	     comd1[3], valtmp, comd1[5], csum1);
   
   snprintf(string, sizeof(string)-1, "%02x %02x %02x %02x %02x %02x %02x",
	    comd1[0], comd1[1], comd1[2], comd1[3], valtmp, comd1[5], csum1);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) command = %s", __FILE__, __LINE__, string);
   
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   /* 
    * Command line 3 
    */
   snprintf(hexstring, sizeof(hexstring)-1, "%08lx",
	    (unsigned long)(40e8/(unsigned long)(count*100)));
   strncpy(strtmp, hexstring + 4, 2);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);
   
   csum2 = comd2[0] ^ comd2[1] ^ comd2[2] ^ comd2[3] ^ valtmp ^ comd2[5];
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) comd= %02x %02x %02x %02x %02x %02x %02x",
	     __FILE__, __LINE__, comd2[0], comd2[1], comd2[2],
	     comd2[3], valtmp, comd2[5], csum2);

   snprintf(string, sizeof(string)-1, "%02x %02x %02x %02x %02x %02x %02x",
	   comd2[0], comd2[1], comd2[2], comd2[3], valtmp, comd2[5], csum2);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) command = %s", __FILE__, __LINE__, string);
   
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   /* 
    * Command line 4 
    */
   snprintf(hexstring, sizeof(hexstring)-1, "%08lx",
	    (unsigned long)(40e8/(unsigned long)(count*100)));
   strncpy(strtmp, hexstring+6, 2);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);
   
   csum3 = comd3[0] ^ comd3[1] ^ comd3[2] ^ comd3[3] ^ valtmp ^ comd3[5];
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) comd= %02x %02x %02x %02x %02x %02x %02x",
	     __FILE__, __LINE__, comd3[0], comd3[1], comd3[2],
	     comd3[3], valtmp, comd3[5], csum3);

   snprintf(string, sizeof(string)-1, "%02x %02x %02x %02x %02x %02x %02x",
	   comd3[0], comd3[1], comd3[2], comd3[3], valtmp, comd3[5], csum3);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) command = %s", __FILE__, __LINE__, string);
   
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   return PASS;
}


/*
 * Get the frame rate for the guide camera 
 */
static PASSFAIL 
getGuiderFrameRate(double *count){
   unsigned long value;
   char hexstring[50]="";
   char longreps[256];
   
   char *string=NULL;
   char *response=NULL;
   char **tokens;
   
   hexstring[0] = '\0';

   
   //snprintf(string, sizeof(string)-1, "53 e0 01 dd 50 3f");
   string="53 e0 01 dd 50 3f";
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   //snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   string="53 e1 01 50 e3";
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   //snprintf(string, sizeof(string)-1, "%s", hexstring);
   //snprintf(hexstring, sizeof(hexstring)-1, "%s%s",
	//    string, *(tokens));
   sprintf(hexstring,"%s%s",hexstring,*(tokens));
   free(tokens);
   
   //snprintf(string, sizeof(string)-1, "53 e0 01 de 50 3c");
   string = "53 e0 01 de 50 3c";
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   //snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   string = "53 e1 01 50 e3";
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   //snprintf(string, sizeof(string)-1, "%s", hexstring);
//    snprintf(hexstring, sizeof(string)-1, "%s%s",
// 	    string, *(tokens));
   sprintf(hexstring,"%s%s",hexstring,*(tokens));
   free(tokens);
   
   //snprintf(string, sizeof(string)-1, "53 e0 01 df 50 3d");
   string = "53 e0 01 df 50 3d";
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   //snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   string = "53 e1 01 50 e3";
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   //snprintf(string, sizeof(string), "%s", hexstring);
//    snprintf(hexstring, sizeof(hexstring)-1, "%s%s",
// 	    string, *(tokens));
   sprintf(hexstring,"%s%s",hexstring,*(tokens));
   free(tokens);
   
   //snprintf(string, sizeof(string)-1, "53 e0 01 e0 50 02");
   string = "53 e0 01 e0 50 02";
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   //snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   string = "53 e1 01 50 e3";
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   //snprintf(string, sizeof(string)-1, "%s", hexstring);
//    snprintf(hexstring, sizeof(hexstring)-1, "%s%s",
// 	    string, *(tokens));
   sprintf(hexstring,"%s%s",hexstring,*(tokens));
   free(tokens);
   
   value = strtoul(hexstring, 0, 16);
   if (value == 0){
      *count = 0;
   } else{
      *count = 40e6/value;
   }

   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) hexstring = %s, count = %lu, Frame rate=%4.2f",
	     __FILE__, __LINE__, hexstring, value, *count);
   
   return PASS;
}


/*
 * Set the exposure time for the camera 
 */
static PASSFAIL
setGuiderExptime(unsigned long count) {

   unsigned char hexstring[20];
   char string[80];
   char *response;
   
   char strtmp[3];
   
   unsigned int comd0[6]={0x53,0xe0,0x02,0xee,0x00,0x50};
   unsigned int comd1[6]={0x53,0xe0,0x02,0xef,0x00,0x50};
   unsigned int comd2[6]={0x53,0xe0,0x02,0xf0,0x00,0x50};
   unsigned int comd3[6]={0x53,0xe0,0x02,0xf1,0x00,0x50};
   
   unsigned int valtmp, csum0, csum1, csum2, csum3;
   
   snprintf(hexstring, sizeof(hexstring)-1, "%08lx", count);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) hexstring = %s count = %lu", 
	     __FILE__, __LINE__, hexstring, count);
   
   /*
    * Preparing serial command 
    */

   /* 
    * Command line 1 
    */
   strncpy(strtmp, hexstring, 2);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);

   csum0 = comd0[0] ^ comd0[1] ^ comd0[2] ^ comd0[3] ^ valtmp ^ comd0[5];
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "comd= %02x %02x %02x %02x %02x %02x %02x",
	     comd0[0], comd0[1], comd0[2], comd0[3], valtmp, comd0[5], csum0);
   
   snprintf(string, sizeof(string)-1, "%02x %02x %02x %02x %02x %02x %02x",
	   comd0[0], comd0[1], comd0[2], comd0[3], valtmp, comd0[5], csum0);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) command = %s", __FILE__, __LINE__, string);
   
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   /* 
    * Command line 2 
    */
   snprintf(hexstring, sizeof(hexstring)-1, "%08lx",
	    count);
   strncpy(strtmp, hexstring + 2, 3);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);

   csum1 = comd1[0] ^ comd1[1] ^ comd1[2] ^ comd1[3] ^ valtmp ^ comd1[5];
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "comd= %02x %02x %02x %02x %02x %02x %02x",
	     comd1[0], comd1[1], comd1[2], comd1[3], valtmp, comd1[5], csum1);
   
   snprintf(string, sizeof(string)-1, "%02x %02x %02x %02x %02x %02x %02x",
	   comd1[0], comd1[1], comd1[2], comd1[3], valtmp, comd1[5], csum1);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) command = %s", __FILE__, __LINE__, string);
   
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   /* 
    * Command line 3 
    */
   snprintf(hexstring, sizeof(hexstring)-1, "%08lx", count);
   strncpy(strtmp, hexstring + 4, 2);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);
   
   csum2 = comd2[0] ^ comd2[1] ^ comd2[2] ^ comd2[3] ^ valtmp ^ comd2[5];
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "comd= %02x %02x %02x %02x %02x %02x %02x",
	     comd2[0], comd2[1], comd2[2], comd2[3], valtmp, comd2[5], csum2);

   snprintf(string, sizeof(string)-1, "%02x %02x %02x %02x %02x %02x %02x",
	    comd2[0], comd2[1], comd2[2], comd2[3], valtmp, comd2[5], csum2);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) command = %s", __FILE__, __LINE__, string);
   
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   /* 
    * Command line 4 
    */
   snprintf(hexstring, sizeof(hexstring)-1, "%08lx", count);
   strncpy (strtmp, hexstring + 6, 2);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);
   
   csum3 = comd3[0] ^ comd3[1] ^ comd3[2] ^ comd3[3] ^ valtmp ^ comd3[5];
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "comd= %02x %02x %02x %02x %02x %02x %02x",
	     comd3[0], comd3[1], comd3[2], comd3[3], valtmp, comd3[5], csum3);

   snprintf(string, sizeof(string)-1, "%02x %02x %02x %02x %02x %02x %02x",
	    comd3[0], comd3[1], comd3[2], comd3[3], valtmp, comd3[5], csum3);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) command = %s", __FILE__, __LINE__, string);
   
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   return PASS;
}


/*
 * Get the exposure time for the camera 
 */
static PASSFAIL 
getGuiderExptime(unsigned long *count) {

   char hexstring[50];
   char longreps[256];
   
   char *response;
   char string[256];
   char** tokens;
   
   hexstring[0] = '\0';
   
   snprintf(string, sizeof(string)-1, "53 e0 01 ee 50 0c");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(string, sizeof(string)-1, "%s", hexstring);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s", string, *(tokens));
   free(tokens);
   
   snprintf(string, sizeof(string)-1, "53 e0 01 ef 50 0d");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(string, sizeof(string)-1, "%s", hexstring);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s", string, *(tokens));
   free(tokens);
   
   snprintf(string, sizeof(string)-1, "53 e0 01 f0 50 12");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(string, sizeof(string)-1, "%s", hexstring);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s", string, *(tokens));
   free(tokens);
   
   snprintf(string, sizeof(string)-1, "53 e0 01 f1 50 13");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(string, sizeof(string)-1, "%s", hexstring);
   snprintf(hexstring, sizeof(hexstring), "%s%s", string, *(tokens));
   free(tokens);
   
   *count = strtoul(hexstring, 0, 16);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) hexstring = %s count = %lu",
	     __FILE__, __LINE__, hexstring, *count);
   
   return PASS;
}


/*
 * Get the digital gain for the camera 
 */
static PASSFAIL
getDigitalGain(int *value) {

   char hexstring[50];
   char longreps[256];
   char string[256];
   
   char *response;
   char **tokens;
   
   hexstring[0] = '\0';
   
   snprintf(string, sizeof(string)-1, "53 e0 01 c6 50 24");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(string, sizeof(string)-1, "%s", hexstring);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s", string, *(tokens));
   free(tokens);
   
   snprintf(string, sizeof(string)-1, "53 e0 01 c7 50 25");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string), "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');

   snprintf(string, sizeof(string)-1, "%s", hexstring);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s", string, *(tokens));
   free(tokens);
   
   *value = strtoul(hexstring, 0, 16) / 256;
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) hexstring = %s gain=%i",
	     __FILE__, __LINE__, hexstring, *value);

   return PASS;
}


/*
 * Set the value of the digital gain
 */
static PASSFAIL 
setDigitalGain(int value){

   unsigned char hexstring[50];
   char string[256];
   char *response;
   
   char strtmp[3];
   
   unsigned int comd0[6]={0x53,0xe0,0x02,0xc6,0x00,0x50};
   unsigned int comd1[6]={0x53,0xe0,0x02,0xc7,0x00,0x50};
   
   unsigned int valtmp, csum0, csum1;
   
   snprintf(hexstring, sizeof(hexstring)-1, "%04lx",
	    (unsigned long)(value*256));
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) hexstring = %s value=%i",
	     __FILE__, __LINE__, hexstring, value);
   
   /*
    * Preparing serial command 
    */

   /* 
    * Command line 1 
    */
   strncpy(strtmp, hexstring, 2);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);
   csum0 = comd0[0] ^ comd0[1] ^ comd0[2] ^ comd0[3] ^ valtmp ^ comd0[5];
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) comd= %02x %02x %02x %02x %02x %02x %02x",
	     __FILE__, __LINE__, comd0[0], comd0[1], comd0[2], comd0[3], 
	     valtmp, comd0[5], csum0);
   
   snprintf(string, sizeof(string)-1, "%02x %02x %02x %02x %02x %02x %02x",
	    comd0[0], comd0[1], comd0[2], comd0[3], valtmp, comd0[5], csum0);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) command=%s", __FILE__, __LINE__, string);
   
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   /* 
    * Command line 2 
    */
   snprintf(hexstring, sizeof(hexstring)-1, "%04lx", 
	    (unsigned long)(value*256));
   strncpy(strtmp, hexstring + 2, 3);
   strtmp[2] = '\0';
   valtmp = strtoul(strtmp, 0, 16);
   
   csum1 = comd1[0] ^ comd1[1] ^ comd1[2] ^ comd1[3] ^ valtmp ^ comd1[5];
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) comd= %02x %02x %02x %02x %02x %02x %02x",
	     __FILE__, __LINE__, comd1[0], comd1[1], comd1[2], comd1[3], 
	     valtmp, comd1[5], csum1);
   
   snprintf(string, sizeof(string)-1, "%02x %02x %02x %02x %02x %02x %02x",
	    comd1[0], comd1[1], comd1[2], comd1[3], valtmp, comd1[5], csum1);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) command=%s", __FILE__, __LINE__, string);
   
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) not able to send command to camera", 
		__FILE__, __LINE__);
      return FAIL;
   }
   
   return PASS;
}


/*
 * Set the guider gain mode
 */
static PASSFAIL 
setGuiderGainMode(int mode){

   char string[256];
   char *response;
   
   if (mode == LOWGAIN){
      snprintf(string, sizeof(string)-1, "53 e0 02 f2 00 50 13");
      if (pdvSerialWriteRead(string, &response) != PASS) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		   "(%s:%d) error occurred sending %s to the camera",
		   __FILE__, __LINE__, string);
	 return FAIL;
      }
      
      if (strcmp ("50 13", response) != 0) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		   "(%s:%d) not able to switch to low gain mode", 
		   __FILE__, __LINE__);
	 return FAIL;
      }
   }
   
   if (mode == HIGHGAIN){
      snprintf(string, sizeof(string), "53 e0 02 f2 06 50 15");
      if (pdvSerialWriteRead(string, &response) != PASS) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		   "(%s:%d) there is no response from the camera", 
		   __FILE__, __LINE__);
	 return FAIL;
      }
      
      if (strcmp ("50 15", response) != 0) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		   "(%s:%d) not able to switch to high gain mode", 
		   __FILE__, __LINE__);
	 return FAIL;
      }
   }
   
   if (mode > 1) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) there is no this specified gain mode for this"
		" camera", __FILE__, __LINE__);
      return FAIL;
   }
   
   return PASS;
}


/* 
 * Get the guider gain mode
 */
static PASSFAIL 
getGuiderGainMode(int *mode) {
   int value;
   char hexstring[50];
   char longreps[256];
   char string[256];
   
   char *response;
   char **tokens;
   
   hexstring[0] = '\0';
   
   snprintf(string, sizeof(string)-1, "53 e0 01 f2 50 10");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(string, sizeof(string)-1, "%s", hexstring);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s", string, *(tokens));
   free(tokens);
   
   value = strtoul(hexstring, 0, 16);
   if (value == 0){
      *mode = LOWGAIN;
   }
   if (value == 6){
      *mode = HIGHGAIN;
   }
   
   return PASS;
}


/* 
 * Check the temperature of the camera
 */
static PASSFAIL
checkGuiderTemp(float *temp) {

   char hexstring[50];
   char longreps[256];
   char string[256];
   
   char *response;
   char **tokens;
   unsigned long adc40d, adc0d;
   
   float value, slope, co;

   hexstring[0] = '\0';
   
   /* 
    * Getting manufacturing data from camera 
    */
   snprintf(string, sizeof(string)-1, "53 ae 05 01 00 00 02 00 50 ab");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 af 12 50 be");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s", 
	    *(tokens+13), *(tokens+12));
   adc40d = strtoul(hexstring, 0, 16);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) ADC40d --> HEX = %s value=%lu",
	     __FILE__, __LINE__, hexstring, adc40d);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s",
	    *(tokens+11), *(tokens+10));
   adc0d = strtoul(hexstring, 0, 16);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) ADC0d --> HEX = %s value=%lu",
	     __FILE__, __LINE__, hexstring, adc0d);
   free(tokens);
   
   /* 
    * Reading Current TEC setpoint 
    */
   snprintf(string, sizeof(string)-1, "53 e0 01 6e 50 8c");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(hexstring, sizeof(hexstring)-1, "%s", *(tokens));
   free(tokens);
   
   snprintf(string, sizeof(string)-1, "53 e0 01 6f 50 8d");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   snprintf(string, sizeof(string)-1, "53 e1 01 50 e3");
   if (pdvSerialWriteRead(string, &response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"(%s:%d) error occurred sending %s to the camera",
		__FILE__, __LINE__, string);
      return FAIL;
   }
   strncpy(longreps, response, sizeof(longreps));
   tokens = stringSplit(longreps, ' ');
   snprintf(string, sizeof(string), "%s", hexstring);
   snprintf(hexstring, sizeof(hexstring)-1, "%s%s",
	    string, *(tokens));
   free(tokens);
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
	     "(%s:%d) HEX = %s", __FILE__, __LINE__, hexstring);
   value = (float)strtoul(hexstring, 0, 16);
   slope = (40.0) / ((float)adc40d - (float)adc0d);
   co = -(slope * (float)adc0d);
   
   *temp = (slope * value) + co;

   return PASS;
}

/*
 * Process the guider image
 */


static PASSFAIL
detrendImage(unsigned char *raw,  unsigned char *flat,  unsigned char *dark,  unsigned char *mask){

	int i,size;
	unsigned char val;


	#ifdef DEBUG
	printf("Begin to processing image....\n");
	#endif

    /* Go through each pixel */
    for (i=0;i<size;i++){
			val=raw[i];
			val=val-dark[i];
			val=val/flat[i];
			val=val/mask[i];
			//if (i == 20400){printf("pixel value=%f\n",val);}

			if (val > SATURATION){
				raw[i]=0;
			}
			if (val < -SATURATION){
				raw[i]=val;
			}

	}

	return PASS;

}


/*
 * Take a pointer to image data and create a FITS image using this data
 * and send it to STDOUT
 */
static PASSFAIL
writeFITSImage(unsigned char *image_p) {

   static unsigned short *image_copy = NULL;
   static int image_width = 0;
   static int image_height = 0;
   HeaderUnit hu;
   time_t date = time(NULL);
   int fd = STDOUT_FILENO;
   char fitscard[FH_MAX_STRLEN];
   char buf[255];
   struct timeval tv;
   struct timezone tz;
   fh_result fh_error;
   int y;


   /*
    * Create the header unit
    */
   hu = fh_create();

   /*
    * Populate the headers 
    */
   fh_set_bool(hu, FH_AUTO, "SIMPLE", 1, "Standard FITS");
   fh_set_int(hu,  FH_AUTO, "BITPIX", 16,"16-bit data");
   fh_set_int(hu,  FH_AUTO, "NAXIS",  2, "Number of axes");
   fh_set_int(hu,  FH_AUTO, "NAXIS1", serv_info->image_width, 
	      "Number of pixel columns");
   fh_set_int(hu,  FH_AUTO, "NAXIS2", serv_info->image_height, 
	      "Number of pixel rows");
   fh_set_int(hu,  FH_AUTO, "PCOUNT", 0, "No 'random' parameters");
   fh_set_int(hu,  FH_AUTO, "GCOUNT", 1, "Only one group");
   
   strftime(fitscard, sizeof(fitscard)-1, "%Y-%m-%dT%T", gmtime(&date));
   fh_set_str(hu, FH_AUTO, "DATE", fitscard, "UTC Date of file creation");
   strftime(fitscard, sizeof(fitscard)-1, "%a %b %d %H:%M:%S %Z %Y", 
	    localtime(&date));
   fh_set_str(hu, FH_AUTO, "HSTTIME", fitscard, "Local time in Hawaii");
   gettimeofday(&tv, &tz);
   fh_set_flt(hu, FH_AUTO, "UNIXTIME", tv.tv_sec + (tv.tv_usec / 1000000.0), 
	      13, "Fractional UNIX timestamp when image was taken");
   fh_set_str(hu, FH_AUTO, "ORIGIN", "CFHT", "Canada-France-Hawaii Telescope");
   fh_set_flt(hu, FH_AUTO, "BZERO", 32768.0, 6,	"Zero factor");
   fh_set_flt(hu, FH_AUTO, "BSCALE", 1.0, 2, "Scale factor");
   fh_set_flt(hu, FH_AUTO, "ETIME", serv_info->exposure_time, 6, 
	      "Integration time (ms)");
   if (serv_info->fits_comment != NULL) {
      fh_set_str(hu, FH_AUTO, "ETYPE", "GUIDE", "Exposure type");
      fh_set_str(hu, FH_AUTO, "IMGINFO", serv_info->fits_comment, 
		 "Sequence details");
   }
   else {
      fh_set_str(hu, FH_AUTO, "ETYPE", "IMAGING", "Exposure type");
      fh_set_str(hu, FH_AUTO, "IMGINFO", "Random imaging sequence", 
		 "Sequence details");
   }
   fh_set_flt(hu, FH_AUTO, "FRMRATE", serv_info->frame_rate, 4,
	      "Requested frame rate (Hz)");
   fh_set_flt(hu, FH_AUTO, "TEMP", serv_info->tec_setpoint, 6, 
	      "TEC cooler setpoint (C)");
   gettimeofday(&tv, &tz);
   fh_set_int(hu, FH_AUTO, "SEQNUM", ++(serv_info->frame_sequence), 
	      "Frame sequence number");
   fh_set_flt(hu, FH_AUTO, "PIXSCALE", PIXSCALE, 5, 
	      "Pixel scale (arcseconds / pixel)");
   if (serv_info->guide_on == TRUE){
	   fh_set_int(hu, FH_AUTO, "WIN_X0", serv_info->win_x0, "X0 coordinate for the camera raster");
	   fh_set_int(hu, FH_AUTO, "WIN_Y0", serv_info->win_y0, "Y0 coordinate for the camera raster");
	   fh_set_int(hu, FH_AUTO, "WIN_X1", serv_info->win_x0 + serv_info->image_width - 1,
	   	      "X1 coordinate for the camera raster");
	      fh_set_int(hu, FH_AUTO, "WIN_Y1", serv_info->win_y0 + serv_info->image_height - 1,
	   	      "Y1 coordinate for the camera raster");
	   fh_set_flt(hu, FH_AUTO, "GD_XOFF", serv_info->guide_xoff,5, "Guide star offset in X");
	   fh_set_flt(hu, FH_AUTO, "GD_YOFF", serv_info->guide_yoff,5, "Guide star offset in Y");

   } else {
	   fh_set_int(hu, FH_AUTO, "WIN_X1", serv_info->image_width - 1,
		      "X1 coordinate for the camera raster");
	   fh_set_int(hu, FH_AUTO, "WIN_Y1", serv_info->image_height - 1,
		      "Y1 coordinate for the camera raster");

   }
   fh_set_flt(hu, FH_AUTO, "NULL_X", (serv_info->image_width - 1) / 2.0, 5, 
	      "Null position (center of aperture hole in X");
   fh_set_flt(hu, FH_AUTO, "NULL_Y", (serv_info->image_height - 1) / 2.0, 5, 
	      "Null position (center of aperture hole in Y");

   /*
    * Extract the telescope position for the image from the Status Server
    * and store this information in the FITS headers.
    */
   //if (ssGetString(SS_TELESCOPE_RA, buf, sizeof(buf) - 1) != PASS) {
   //   cfht_logv(CFHT_MAIN, CFHT_DEBUG,
	//	"(%s:%d) unable to retrieve %s from the"
	//	" Status Server: %s\n", __FILE__, __LINE__, 
	//	SS_TELESCOPE_RA, ssGetStrError());
    //  fh_set_str(hu, FH_AUTO, "RA", "UNKNOWN", "Telescope RA");
   //}
   //else {
      fh_set_str(hu, FH_AUTO, "RA", buf, "Telescope RA");
   //}
   //if (ssGetString(SS_TELESCOPE_DEC, buf, sizeof(buf) - 1) != PASS) {
   //   cfht_logv(CFHT_MAIN, CFHT_DEBUG, 
//		"(%s:%d) unable to retrieve %s from the"
//		" Status Server: %s\n", __FILE__, __LINE__, 
//		SS_TELESCOPE_DEC, ssGetStrError());
//      fh_set_str(hu, FH_AUTO, "DEC", "UNKNOWN", "Telescope Dec");
//   }
//   else {
      fh_set_str(hu, FH_AUTO, "DEC", buf, "Telescope Dec");
//   }
//   if (ssGetString(SS_TELESCOPE_EQ, buf, sizeof(buf) - 1) != PASS) {
//      cfht_logv(CFHT_MAIN, CFHT_DEBUG, 
//		"(%s:%d) unable to retrieve %s from the"
//		" Status Server: %s\n", __FILE__, __LINE__, 
//		SS_TELESCOPE_EQ, ssGetStrError());
      fh_set_str(hu, FH_AUTO, "EQUINOX", "UNKNOWN", "Telescope Equinox");
//   }
//   else {
      fh_set_str(hu, FH_AUTO, "EQUINOX", buf, "Telescope Equinox");
//   }

   /*
    * If the frame count has been reached, clear out the save 
    * information
    */
   if ((serv_info->fits_comment != NULL) &&
       (serv_info->frame_sequence >= serv_info->frame_save_count)) {
      free(serv_info->fits_comment);
      serv_info->fits_comment = NULL;
      serv_info->frame_save_count = 0;
      serv_info->frame_sequence = 0;
   }

   /* 
    * Write out the FITS header 
    */
   if ((fh_error = fh_write(hu, fd)) != FH_SUCCESS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) unable to write FITS header"
		" (fh_error = %d)\n", __FILE__, __LINE__, fh_error);
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY, 
		"%s (errno=%d)", strerror(errno), errno);
      return FAIL;
   }

   /*
    * Allocate memory to store the final (correctly flipped) image
    * if it hasn't been allocated yet.
    */
   if (serv_info->image_width != image_width ||
       serv_info->image_height != image_height) {
      if (image_copy != NULL) {
	 free(image_copy);
      }
      image_width = serv_info->image_width;
      image_height = serv_info->image_height;

      image_copy = (unsigned short *)cli_malloc(image_width * 
						image_height *
 						sizeof(uint16_t));
   }

   /*
    * Copy and flip the image vertically
    */
   for (y = 0; y < image_height; y++) {
      memcpy(&(image_copy[y * image_width]), 
	     &(((unsigned short *)(image_p))[(image_height - y - 1) *
					     image_width]),
	       image_width * sizeof(uint16_t));
   }

   /*
    * Write out the image data
    */
//   if ((fh_error = fh_write_padded_image(hu, fd, (unsigned short *)image_p, 
//					 serv_info->image_width *
//					 serv_info->image_height *
//					 sizeof(uint16_t), FH_TYPESIZE_16U)) 
//       != FH_SUCCESS) {
   if ((fh_error = fh_write_padded_image(hu, fd, 
					 (unsigned short *)image_copy, 
					 image_width * image_height *
					 sizeof(uint16_t), FH_TYPESIZE_16U)) 
       != FH_SUCCESS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) unable to write FITS image data"
		" (fh_error = %d)\n", __FILE__, __LINE__, fh_error);
      return FAIL;
   }

   /*
    * Free up the memory for the FITS header
    */
   fh_destroy(hu);

   return PASS;
}


/*
 * Handle a new client connection
 */
static void *
clientAdd(unsigned char remote_ip[4])
{
   client_info_t *cinfo;
   struct hostent *hp;	  /* Host entry returned from gethostbyaddr */

   /*
    * Allocate memory for the client data structure
    */
   cinfo = (client_info_t *)cli_malloc(sizeof(*cinfo));
   memset(cinfo, 0, sizeof(client_info_t));
   memcpy(cinfo->remote_ip, remote_ip, sizeof (cinfo->remote_ip));

   /* 
    * Determine the hostname from the IP address 
    */
   hp = gethostbyaddr(cinfo->remote_ip, sizeof (struct in_addr), AF_INET);
   if (hp == NULL) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) gethostbyaddr call failed", __FILE__, __LINE__);
      cinfo->hostname = cli_strdup("UNKNOWN");
   }
   else {
      cinfo->hostname = cli_strdup(hp->h_name);
   }

   /* 
    * Store the time the client connected to the Control Server.  This
    * is logged when a request is made to list the connected clients.
    */
   time(&(cinfo->connect_ts));

   /* 
    * Add the client to the linked list of clients connected to the 
    * server.
    */
   appendDataToList(cinfo, serv_info->client_list, cli_malloc);

   return cinfo;
}


/*
 * Handle a client disconnect
 */
static void
clientDelete(void *client, char *buffer)
{
   client_info_t *cinfo = (client_info_t *)client;

   /* 
    * Search the client linked list and remove the client object
    * from the list.
    */
   searchListDeleteData(cinfo, serv_info->client_list, comparePointer);

   /*
    * Free up the memory associated with the client object.
    */
   deallocateClient(cinfo);
}


/*
 * Handle a command request from a client
 */
static void
clientReceive(void *client, char *buffer)
{
   int cargc = 0;               /* number of arguments found in args */
   static char **cargv = NULL;  /* arguments in split form */
//   client_info_t *cinfo = (client_info_t *)client;
   char *buf_p;
   char *p;

   /* This block is related to the ISU management */
   BOOLEAN homed;

   /*
    * Log the message received by the server
    */
   cfht_logv(CFHT_MAIN, CFHT_DEBUG,
	     "(%s:%d) RECV> %s", __FILE__, __LINE__, buffer);
   
   /* 
    * Advance past the command for argument parsing purposes.
    */
   buf_p = trim(buffer);
   p = strchr(buf_p, ' ');
   if (p == NULL) {

      /* ----------------------------------------------------
       * Handle all commands which don't have any parameters
       * ----------------------------------------------------*/

      /*
       * Handle the various ways the client can disconnect from the 
       * Control Server.
       */
      if (!strcasecmp(buf_p, QUIT_CMD) ||
	  !strcasecmp(buf_p, BYE_CMD) ||
	  !strcasecmp(buf_p, EXIT_CMD) || 
	  !strcasecmp(buf_p, LOGOUT_CMD)) {
	 buffer[0] = '\0';
	 return;
      }

      /*
       * Shutdown the Raptor guide camera server. 
       */
      if (!strcasecmp(buf_p, SHUTDOWN_CMD)) {
	 buffer[0] = '\0';
	 serv_info->serv_done = 1;
	 return;
      }

      /*
       * Handle a query of the current exposure time
       */
      if (!strcasecmp(buf_p, EXPTIME_CMD)) {
	 unsigned long count;
	 
	 if (getGuiderExptime(&count) != PASS) {
	    cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		      "(%s:%d) unable to read exposure time from the"
		      " Raptor camera", __FILE__, __LINE__);
	    sprintf(buffer, "%c %s \"Unable to exposure time in the camera\"", 
		    FAIL_CHAR, EXPTIME_CMD);
	    cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		      "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	    return;
	 }
	 serv_info->exposure_time = count / 40e3;

	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) exposure time is %.9f ms",
		   __FILE__, __LINE__, serv_info->exposure_time);
	 sprintf(buffer, "%c %s %.9f", PASS_CHAR, EXPTIME_CMD, 
		 serv_info->exposure_time);
	 
	 return;
      }

      /*
       * Handle a query of the frame rate
       */
      if (!strcasecmp(buf_p, FRAMERATE_CMD)) {
    	  double frame_rate;
    	  char* camera_response=NULL;
    	  sprintf(buffer, "%c %s", PASS_CHAR, FRAMERATE_CMD);
    	  cfht_logv(CFHT_MAIN, CFHT_DEBUG,
 		      "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);

      /* Make sure the connection with camera is still alive */
	  if (checkCameraStatus(&camera_response) != PASS) {
		  cfht_logv(CFHT_MAIN, CFHT_ERROR,
			"(%s:%d) there is no response from the camera when checking"
			" the camera status - exiting", __FILE__, __LINE__);
	  }

 	 if (getGuiderFrameRate(&frame_rate) != PASS) {
    	    cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
 		      "(%s:%d) unable to read frame rate from the"
 		      " Raptor camera", __FILE__, __LINE__);
 	     sprintf(buffer, "%c %s \"Unable to read frame rate in"
 		    " the camera\"", FAIL_CHAR, FRAMERATE_CMD);
 	    cfht_logv(CFHT_MAIN, CFHT_DEBUG,
 		      "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
        return;
 	 
 	 }
 	 serv_info->frame_rate = frame_rate;
 	 
 	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
 		   "(%s:%d) frame rate in camera set to %5.2f Hz",
 		   __FILE__, __LINE__, frame_rate);
	 sprintf(buffer, "%c %s %5.2f", PASS_CHAR, FRAMERATE_CMD, frame_rate);
	 
	 return;
      }

      /*
       * Handle a query for the TEC set point
       */
      if (!strcasecmp(buf_p, TEC_CMD)) {
	 float setpoint;

	 if (getGuiderTECPoint(&setpoint) != PASS) {
	    cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		      "(%s:%d) unable to read TEC set point from the"
		      " Raptor camera", __FILE__, __LINE__);
	    sprintf(buffer, "%c %s \"Unable to read TEC set point in"
		    " the camera\"", FAIL_CHAR, TEC_CMD);
	    cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		      "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	    return;
	 }
	 serv_info->tec_setpoint = setpoint;
	 
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) TEC set point is %4.1f degrees C",
		   __FILE__, __LINE__, setpoint);
	 sprintf(buffer, "%c %s %4.1f degrees C", PASS_CHAR, 
		 TEC_CMD, setpoint);
	 
	 return;
      }

      /*
       * Handle a query for the temperature
       */
      if (!strcasecmp(buf_p, TEMP_CMD)) {
	 float temp;

	 if (checkGuiderTemp(&temp) != PASS) {
	    cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		      "(%s:%d) unable to read temperature from the"
		      " Raptor camera", __FILE__, __LINE__);
	    sprintf(buffer, "%c %s \"Unable to read temperature in"
		    " the camera\"", FAIL_CHAR, TEMP_CMD);
	    cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		      "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	    return;
	 }
	 serv_info->temp = temp;
	 
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) Temperarure is %4.1f degrees C",
		   __FILE__, __LINE__, temp);
	 sprintf(buffer, "%c %s %4.1f degrees C", PASS_CHAR, 
		 TEMP_CMD, temp);
	 
	 return;
      }

      /*
       * Handle a query for the image ROI
       */
      if (!strcasecmp(buf_p, ROI_CMD)) {

          cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
              "(%s:%d) Image ROI is %i X %i now.\n",
    		  __FILE__, __LINE__, pdv_get_width(serv_info->pdv_p),
			  pdv_get_height(serv_info->pdv_p));
    	  sprintf(buffer, "%c %s is %i X %i ", PASS_CHAR,
    		  ROI_CMD, pdv_get_width(serv_info->pdv_p),
			  pdv_get_height(serv_info->pdv_p));

    	  return;
      }
      /*
       * If we made it this far, this is an unrecognized command request
       * from the client which doesn't have parameters.
       */
      sprintf(buffer, "%c \"Syntax Error\"", FAIL_CHAR);
      cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		"(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
      return;
   }

   /* ---------------------------------------------
    * Handle all commands which contain parameters
    * ---------------------------------------------*/
   *p = '\0';

   /*
    * Free up memory associated with arguments previously parsed 
    */
   if (cargv != NULL) {
      cli_argv_free(cargv);
   }
   cargv = cli_argv_quoted(&cargc, ++p);

   /*
    * Handle a frame rate command from a client 
    */
   if (!strcasecmp(buf_p, FRAMERATE_CMD)) {
      double frame_rate;
      char *camera_response = NULL;
      char *stop_at = NULL;

      frame_rate = strtod(cargv[0], &stop_at);
      
      /* Make sure the frame rate is valid */
      if ((errno == ERANGE) || (errno == EINVAL) || (*stop_at != '\0')) {
	 sprintf(buffer, "%c %s \"Invalid Argument Specified\"", 
		 FAIL_CHAR, FRAMERATE_CMD);
	 cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		   "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	 return;
      }
      if (frame_rate <= 0 || (1e3 / frame_rate) > USER_TIMEOUT ||
	  frame_rate > 120) {
	 sprintf(buffer, "%c %s \"Frame Rate Specified is Invalid\"", 
		 FAIL_CHAR, FRAMERATE_CMD);
	 cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		   "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	 return;
      }
      
      /* Make sure the connection with camera is still alive */
      if (checkCameraStatus(&camera_response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_ERROR,
		"(%s:%d) there is no response from the camera when checking"
		" the camera status - exiting", __FILE__, __LINE__);
      }
      /* 
       * Try to update the frame rate of the camera
       */
      //cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		//"(%s:%d) attempting to set the frame rate to %5.2f Hz",
		//__FILE__, __LINE__, frame_rate);
      if (setGuiderFrameRate(frame_rate) != PASS) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) unable to set frame rate to %5.2f Hz",
		   __FILE__, __LINE__, frame_rate);
	 sprintf(buffer, "%c %s \"Unable to set frame rate in the camera\"", 
		 FAIL_CHAR, FRAMERATE_CMD);
	 cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		   "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	 return;
      }

      /*
       * Just to be sure read back the frame rate from the camera
       */
       if (getGuiderFrameRate(&frame_rate) != PASS) {
 	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
 		   "(%s:%d) unable to read frame rate from the"
 		   " Raptor camera", __FILE__, __LINE__);
 	 sprintf(buffer, "%c %s \"Unable to read frame rate in the camera\"", 
 		 FAIL_CHAR, FRAMERATE_CMD);
 	 cfht_logv(CFHT_MAIN, CFHT_DEBUG,
 		   "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
 	 return;
       }
      serv_info->frame_rate = frame_rate;

      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) frame rate in camera set to %5.2f Hz",
		__FILE__, __LINE__, frame_rate);
      sprintf(buffer, "%c %s %5.2f", PASS_CHAR, FRAMERATE_CMD, frame_rate);

      return;
   }

   /*
    * Handle a exposure time command from a client 
    */
   if (!strcasecmp(buf_p, EXPTIME_CMD)) {
      double exptime;
      char *stop_at = NULL;

      exptime = strtod(cargv[0], &stop_at);
      
      /* Make sure the exposure time is valid */
      if ((errno == ERANGE) || (errno == EINVAL) || 
	  (exptime <= 0) || (*stop_at != '\0')) {
	 sprintf(buffer, "%c %s \"Invalid Argument Specified\"", 
		 FAIL_CHAR, EXPTIME_CMD);
	 cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		   "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	 return;
      }
      
      /* 
       * Try to update the exposure time of the camera
       */
      exptime *= 40e3;
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) attempting to set the exposure time to %5.2f ms",
		__FILE__, __LINE__, exptime);
      if (setGuiderExptime(exptime) != PASS) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) unable to set exposure time to %5.2f ms",
		   __FILE__, __LINE__, exptime / 40e3);
	 sprintf(buffer, "%c %s \"Unable to set exposure time in the"
		 " camera\"", FAIL_CHAR, EXPTIME_CMD);
	 cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		   "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	 return;
      }
      serv_info->exposure_time = exptime / 40e3;

      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) exposure time in camera set to %.3f ms",
		__FILE__, __LINE__, serv_info->exposure_time);
      sprintf(buffer, "%c %s %5.2f", PASS_CHAR, EXPTIME_CMD,
	      serv_info->exposure_time);

      return;
   }

   /*
    * Handle a TEC cooler setpoint command 
    */
  if (!strcasecmp(buf_p, TEC_CMD)) {
      float setpoint;
      char *stop_at = NULL;

      setpoint = strtod(cargv[0], &stop_at);
      
      /* Make sure the set point is valid */
      if ((errno == ERANGE) || (errno == EINVAL) || (*stop_at != '\0')) {
	 sprintf(buffer, "%c %s \"Invalid Argument Specified\"", 
		 FAIL_CHAR, TEC_CMD);
	 cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		   "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	 return;
      }
      
      /* 
       * Try to update the TEC set point of the camera
       */
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) attempting to set the TEC setpoint to %5.2f degrees"
		" C", __FILE__, __LINE__, setpoint);
      if (setGuiderTECPoint(setpoint) != PASS) {
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		   "(%s:%d) unable to set TEC setpoint to %5.2f degrees C",
		   __FILE__, __LINE__, setpoint);
	 sprintf(buffer, "%c %s \"Unable to set TEC set point in the"
		 " camera\"", FAIL_CHAR, TEC_CMD);
	 cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		   "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	 return;
      }

      serv_info->tec_setpoint = setpoint;

      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) TEC cooler set point in camera set to %5.2f"
		" degrees C", __FILE__, __LINE__, setpoint);
      sprintf(buffer, "%c %s %5.2f", PASS_CHAR, TEC_CMD, setpoint);

      return;
   }

   /*
    * Handle a request to turn on or off the video from a client 
    */
   if (!strcasecmp(buf_p, VIDEO_CMD)) {
      if ((!strcasecmp(cargv[0], "ON")) && (cargc == 1)) {
	 serv_info->video_on = TRUE;
	 sprintf(buffer, "%c ON", PASS_CHAR);
      } else if ((!strcasecmp(cargv[0], "OFF")) && (cargc == 1)) {
	 serv_info->video_on = FALSE;
	 sprintf(buffer, "%c OFF", PASS_CHAR);
      }
      else {
	 sprintf(buffer, "%c \"Invalid video mode\"", FAIL_CHAR);
      }
      cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		"(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
      return;
   }

   /*
    * Handle a request to turn on or off the ISU mode from a client
    */
   if (!strcasecmp(buf_p, ISU_CMD)) {
      if ((!strcasecmp(cargv[0], "ON")) && (cargc == 1)) {
         /* Check if the ISU needs to perform a new homing */
         if (check_isu_homed(&homed) != PASS) {
	     cfht_logv(CFHT_MAIN, CFHT_WARN, "(%s:%d) Failed checking "
             "homing status of the ISU mechanism",
	      __FILE__, __LINE__);
             homed = 0;
         }  
         if (homed == 0) {
            /* Home the ISU */
            if (home_isu() != PASS) {
	        cfht_logv(CFHT_MAIN, CFHT_ERROR, "(%s:%d) Fatal error on "
                "the fast guiding loop: failed homing the ISU mechanism",
	         __FILE__, __LINE__);
                 exit(EXIT_FAILURE);
            }  
         }
         serv_info->isu_on = TRUE;
	 sprintf(buffer, "%c ON", PASS_CHAR);
#ifdef SIM_STAR
	 cfht_logv(CFHT_MAIN, CFHT_LOGONLY, "(%s:%d) Beware that "
         "SIM_STAR is defined and the centroid algorithm is NOT called",
	 __FILE__, __LINE__);
#endif
#ifdef DEBUG
       cfht_logv(CFHT_MAIN, CFHT_LOGONLY, "(%s:%d) Beware that "
         "DEBUG is defined and status / setup positions are written in "
         "Outputs.csv file", __FILE__, __LINE__);
#endif

      } else if ((!strcasecmp(cargv[0], "OFF")) && (cargc == 1)) {
         /* Stop the ISU */
         if (stop_isu() != PASS) {
	     cfht_logv(CFHT_MAIN, CFHT_WARN, "(%s:%d) Failed stopping "
             " the ISU mechanism", __FILE__, __LINE__);
         }
    	 serv_info->isu_on = FALSE;
    	 sprintf(buffer, "%c OFF", PASS_CHAR);
      }
      else {
    	   sprintf(buffer, "%c \"Invalid guide mode\"", FAIL_CHAR);
      }
      cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		"(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
      return;
   }

   /*
    * Handle a request to turn on or off the guide mode from a client
    */
   if (!strcasecmp(buf_p, GUIDE_CMD)) {
	   if ((!strcasecmp(cargv[0], "ON")) && (cargc == 1)) {
		   serv_info->guide_on = TRUE;
		   sprintf(buffer, "%c ON", PASS_CHAR);
      } else if ((!strcasecmp(cargv[0], "OFF")) && (cargc == 1)) {
    	   serv_info->guide_on = FALSE;
    	   sprintf(buffer, "%c OFF", PASS_CHAR);
      }
      else {
    	   sprintf(buffer, "%c \"Invalid guide mode\"", FAIL_CHAR);
      }
      cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		"(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
      return;
   }

   /*
    * Handle a request to save a sequence of images.  This will set the 
    * FITS headers in such a way that a downstream program saving the data
    * can key on this and save the proper information
    */
   if (!strcasecmp(buf_p, SAVE_CMD)) {
      int save_count;
      char *stop_at = NULL;

      save_count = strtol(cargv[0], &stop_at, 10);
      
      /* Make sure the save cound is valid */
      if ((errno == ERANGE) || (errno == EINVAL) || (*stop_at != '\0') ||
	  (save_count < 1 || save_count > MAX_SAVE_COUNT)) {
	 sprintf(buffer, "%c %s \"Invalid Argument Specified\"", 
		 FAIL_CHAR, SAVE_CMD);
	 cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		   "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	 return;
      }
      serv_info->fits_comment = cli_strdup(cargv[1]);
      serv_info->frame_save_count = save_count;
      serv_info->frame_sequence = 0;
      cfht_logv(CFHT_MAIN, CFHT_DEBUG, 
		"(%s:%d) save count = %d, comment = %s",
		__FILE__, __LINE__, save_count, cargv[1]);
      sprintf(buffer, "%c %s", PASS_CHAR, SAVE_CMD);
      cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		"(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
      return;
   }

   /*
    * Handle a request to set a ROI of images.  This will set the
    * PDV card to change the ROI.
    */
   if (!strcasecmp(buf_p, ROI_CMD)) {
	   int x0,y0,x1,y1;
	   int s_width,s_height;
	   char *stop_at = NULL;

	   if (cargc != 5 ){
		   sprintf(buffer, "%c %s \"Input FOUR numbers to set ROI.\"",
			  FAIL_CHAR, ROI_CMD);
		   cfht_logv(CFHT_MAIN, CFHT_DEBUG,
			  "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
	   }

	   x0 = strtol(cargv[0], &stop_at,10);
	   y0 = strtol(cargv[1], &stop_at,10);
	   x1 = strtol(cargv[2], &stop_at,10);
	   y1 = strtol(cargv[3], &stop_at,10);

	   /* Make sure the set point is valid */
	   if ((errno == ERANGE) || (errno == EINVAL) || (*stop_at != '\0')) {
		   sprintf(buffer, "%c %s \"Invalid Argument Specified\"",
			  FAIL_CHAR, ROI_CMD);
		   cfht_logv(CFHT_MAIN, CFHT_DEBUG,
			  "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
		 return;
	   }

	   /* Adding serv_info->fits_comment, so that this variable won't NULL.  This will enable 
	    *   FITS keyword ETYPE=GUIDE in line 2300 
	    */
	   serv_info->fits_comment = cli_strdup("GUIDING");

	   s_width = (x1 - x0)+1;
	   s_height = (y1 - y0)+1;

	   if (s_height<2 && s_width<2){
		   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
				 "(%s:%d) image size incorrect, check ROI setting",
				   __FILE__, __LINE__);
		   sprintf(buffer, "%c %s \"image size incorrect, check ROI setting",
				 FAIL_CHAR, ROI_CMD);
		   cfht_logv(CFHT_MAIN, CFHT_DEBUG,
				 "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
		   return;
	   }
	   if (pdv_set_roi(serv_info->pdv_p,  x0, s_width,y0, s_height) != 0){
		   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
				 "(%s:%d) unable to set image ROI",
				   __FILE__, __LINE__);
		   sprintf(buffer, "%c %s \"unable to set image ROI\"",
				 FAIL_CHAR, ROI_CMD);
		   cfht_logv(CFHT_MAIN, CFHT_DEBUG,
				 "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
		   return;

	   }
	   if (pdv_enable_roi(serv_info->pdv_p, 1) != 0){
		   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
				 "(%s:%d) set ROI failed",
				   __FILE__, __LINE__);
		   sprintf(buffer, "%c %s \"set ROI failed\"",
				 FAIL_CHAR, ROI_CMD);
		   cfht_logv(CFHT_MAIN, CFHT_DEBUG,
				 "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
		   return;

		}
	   serv_info->image_width=s_width;
	   serv_info->image_height=s_height;
	   serv_info->win_x0=x0;
	   serv_info->win_y0=y0;

	   cfht_logv(CFHT_MAIN, CFHT_DEBUG,
		   "(%s:%d) image ROI setting is: x0=%i y0=%i x1=%i y1=%i",
		   __FILE__, __LINE__, x0, y0, x1, y1);
	   sprintf(buffer, "%c %s", PASS_CHAR, ROI_CMD);
	   cfht_logv(CFHT_MAIN, CFHT_DEBUG,
			"(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);

	   return;
   }
   /*
    * If we made it this far, this is an unrecognized command request
    * from the client.
    */
   sprintf(buffer, "%c \"Syntax Error\"", FAIL_CHAR);
   cfht_logv(CFHT_MAIN, CFHT_DEBUG,
	     "(%s:%d) SEND> %s", __FILE__, __LINE__, buffer);
}


/*
 * Handle a cleanup of the socket resources and make sure the shutter is 
 * closed.
 */
static void
cleanup(void) {

   /*
    * Close the connection to the camera link board
    */
//   if (serv_info->pdv_p != NULL) {
//      pdv_close(serv_info->pdv_p);
//   }

   /*
    * Cleanup the listening socket 
    */
   if (serv_info->raptor_serv != NULL) {
      sockserv_destroy(serv_info->raptor_serv);
      serv_info->raptor_serv = NULL;
   }

   exit(EXIT_SUCCESS);
}


int
main(int argc, const char* argv[])
{
   int edt_debug_level;
   int edt_unit = 0;
   int edt_channel = 0;
   char edt_devname[256];
   char *edt_unitstr = "0";
   Edtinfo edtinfo;
   char bitdir[256];
   char *camera_response;
   int gain_mode=1;
   int digital_gain;
   BOOLEAN last_video_on_state = FALSE;
   unsigned char *image_p;
   int last_timeouts = 0;
   int timeouts;

   /* This block is related to the ISU management */
   float xc = 0;
   float yc = 0;
   BOOLEAN x_fault, y_fault;
   float xangle = 0;
   float yangle = 0;
   double last_x_angle, last_y_angle;
   double next_x_angle, next_y_angle;
   thdata thread_data;
   pthread_t set_analog_slope_t;

#ifdef DEBUG
   double time_spent;
   struct timeval last_time, this_time;
   struct timezone tz;
   int index = 0;
   static FILE* filePtr = NULL;
   double meas_setup_x_angle, meas_setup_y_angle;

   gettimeofday(&last_time,&tz);
#endif

   /*
    * Set up the environment variable used by the cfht_log system to 
    * determine whether debug messages will be logged.
    */
   setenv("CFHTDEBUG", "On", 1);
   
   /*
    * Set up all memory allocations to retry if memory is available.  This 
    * has the added advantage that there is no need to change for memory
    * allocation failures.  If memory isn't available, it will hang until
    * memory becomes available.
    */
   cli_malloc_retry(TRUE);

   /* 
    * Allocate and initialize the server information structure 
    */
   serv_info = (server_info_t *)cli_malloc(sizeof(server_info_t));
   memset(serv_info, 0, sizeof(server_info_t));

   /*
    * Initialize the CFHT logging stuff.
    */
   cfht_log(CFHT_MAIN, CFHT_LOG_ID, argv[0]);
   cfht_logv(CFHT_MAIN, CFHT_START, "%s", argv[0]);

   /*
    * Log on to the Status Server
    */
//   if (ssLogon(argv[0]) != PASS) {
//      cfht_logv(CFHT_MAIN, CFHT_ERROR,/
//		"(%s:%d) Status Server logon failed: %s",
//		__FILE__, __LINE__, ssGetStrError());
      //exit(EXIT_FAILURE);
 //  }

   /*
    * TODO: Look at adding a function callback for logging
    */

   /*
    * Initialize the verbosity level for messages from the Camera Link library
    */
   edt_debug_level = edt_msg_default_level();
   edt_debug_level |= EDTAPP_MSG_INFO_1;
   edt_debug_level |= PDVLIB_MSG_INFO_1;
   edt_debug_level |= PDVLIB_MSG_WARNING;
   edt_debug_level |= PDVLIB_MSG_FATAL;

   /*
    * Additional debug level messages
    */
//   edt_debug_level |= EDTAPP_MSG_INFO_2;
//   edt_debug_level |= PDVLIB_MSG_INFO_2;

   /*
    * Apply the debug level
    */
   edt_msg_set_level(edt_msg_default_handle(), edt_debug_level);

   /*
    * If porting this code to an application, be sure to free this and
    * reallocate if you call pdv_initcam multiple times.
    */
   if ((serv_info->dd_p = pdv_alloc_dependent()) == NULL) {
      cfht_logv(CFHT_MAIN, CFHT_ERROR,
		"(%s:%d) pdv_alloc_dependent() failed - exiting",
		__FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

   /*
    * Read the camera config file
    */
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) using camera config file = %s",
	     __FILE__, __LINE__, RAPTOR_CONFIG);
   if (pdv_readcfg(RAPTOR_CONFIG, serv_info->dd_p, &edtinfo) != 0){
      cfht_logv(CFHT_MAIN, CFHT_ERROR,
		"(%s:%d) pdv_readcfg() failed - exiting",
		__FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

   /*
    * Open the device
    */
   edt_unit = edt_parse_unit_channel(edt_unitstr, edt_devname, 
				     "pdv", &edt_channel);
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) opening %s unit %d", __FILE__, __LINE__,
	     edt_devname, edt_unit);
   if ((serv_info->edt_p = edt_open_channel(edt_devname, 
					    edt_unit, edt_channel)) == NULL) {
      cfht_logv(CFHT_MAIN, CFHT_ERROR,
		"(%s:%d) edt_open(%s%d) failed - exiting",
		__FILE__, __LINE__, edt_devname, edt_unit);
      exit(EXIT_FAILURE);
   }

   /*
    * Initialize the framegrabber board and camera
    */
   if (pdv_initcam(serv_info->edt_p, serv_info->dd_p, edt_unit, 
		   &edtinfo, RAPTOR_CONFIG, bitdir, 0) != 0){
      cfht_logv(CFHT_MAIN, CFHT_ERROR,
		"(%s:%d) pdv_initcam() failed - exiting",
		__FILE__, __LINE__);
      edt_close(serv_info->edt_p);
      exit(EXIT_FAILURE);
   }

   /*
    * Set the clock frequency in MHz on the board
    */
   pdv_cls_set_clock(serv_info->edt_p, 40.0);

   /*
    * Set a reasonable image timeout value based on the image size,
    * exposure time (if set) and pixel clock speed (if set)
    */
   if (pdv_auto_set_timeout(serv_info->edt_p) != 0){
      edt_close(serv_info->edt_p);
      cfht_logv(CFHT_MAIN, CFHT_ERROR,
		"(%s:%d) can not set the camera time out - exiting",
		__FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

   /*
    * TODO: The following statement is part of the guiderinit.c program.
    * I'm not sure it makes sense for the interface to be closed, but I'll
    * leave it here until the testing phase.
    */
   edt_close(serv_info->edt_p);
   cfht_logv(CFHT_MAIN, CFHT_ERROR,
	     "(%s:%d) edt_close() performed", __FILE__, __LINE__);

   /*
    * Check the system status from the camera
    */
   if (checkCameraStatus(&camera_response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_ERROR,
		"(%s:%d) there is no response from the camera when checking"
		" the camera status - exiting", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) camera status OK", __FILE__, __LINE__);

   /* 
    * Turn off non-uniform correction
    */
   if (setGuiderNUC(0) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_ERROR,
		"(%s:%d) not able to set NUC status - exiting", 
		__FILE__, __LINE__);
      exit(EXIT_FAILURE);
      
   }
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) camera NUC turned off", __FILE__, __LINE__);

   /* 
    * Turn off Auto Level 
    */
   if (setGuiderAutoLevel(0) != 0) {
      cfht_logv(CFHT_MAIN, CFHT_ERROR,
		"(%s:%d) not able to set camera auto level - exiting", 
		__FILE__, __LINE__);
      exit(EXIT_FAILURE);
      
   }
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) camera auto level turned off", __FILE__, __LINE__);

   /*
    * Enable TEC cooler
    */
   if (enableGuiderTEC() != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) enable TEC control failed - exiting", 
		__FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) camera TEC enabled", 
	     __FILE__, __LINE__);

   if (checkCameraStatus(&camera_response) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_ERROR,
		"(%s:%d) there is no response from the camera when checking"
		" the camera status - exiting", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

   /*
    * Set the camera to the high gain mode
    */
   if (setGuiderGainMode(HIGHGAIN) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) unable to set camera to high-gain mode - exiting",
		__FILE__, __LINE__);
      //exit(EXIT_FAILURE);
   }

   /*
    * Read back the gain mode to verify that the guider is in the high
    * gain mode
    */
   if (getGuiderGainMode(&gain_mode) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) unable to retrieve gain mode from the"
		" camera - exiting", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   if (gain_mode != HIGHGAIN) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) gain read back from the camera=%d which is not"
		" the expected high-gain mode - exiting",
		__FILE__, __LINE__, gain_mode);
      exit(EXIT_FAILURE);
   }
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) camera is set to HIGH gain mode",
	     __FILE__, __LINE__);

   /*
    * Set a temperature set point
    */
   if (setGuiderTECPoint(DEFAULT_TEC_SETPOINT) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) unable to set the TEC setpoint to %.2f"
		" degrees - exiting", __FILE__, __LINE__, 
		DEFAULT_TEC_SETPOINT);
      exit(EXIT_FAILURE);
   }
   serv_info->tec_setpoint = DEFAULT_TEC_SETPOINT;
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) TEC setpoint set to %.2f degrees",
	     __FILE__, __LINE__, DEFAULT_TEC_SETPOINT);

   
   /*
    * Set the digital gain value to be 1
    */
   if (setDigitalGain(DEFAULT_DIGITAL_GAIN) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) unable to set default digital gain of %d - exiting",
		__FILE__, __LINE__, DEFAULT_DIGITAL_GAIN);
      exit(EXIT_FAILURE);
   }

   /*
    * Read back the digital gain value and make sure it matches what is 
    * expected
    */
   if (getDigitalGain(&digital_gain) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) unable to retrieve digital gain from the"
		" camera - exiting", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   if (digital_gain != DEFAULT_DIGITAL_GAIN) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) digital gain retrieved from the camera does not"
		" match what is expected (%d != %d) - exiting", 
		__FILE__, __LINE__, digital_gain, DEFAULT_DIGITAL_GAIN);
      exit(EXIT_FAILURE);
   }

   /*
    * Set the default exposure time to be 10 milliseconds
    */
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) attempting to set the exposure time to %.3f ms",
	     __FILE__, __LINE__, DEFAULT_EXPOSURE_TIME);
   if (setGuiderExptime(DEFAULT_EXPOSURE_TIME * 40e3) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) unable to set exposure time to %.3f ms - exiting",
		__FILE__, __LINE__, DEFAULT_EXPOSURE_TIME);
      exit(EXIT_FAILURE);
   }
   serv_info->exposure_time = DEFAULT_EXPOSURE_TIME;
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) exposure time in camera set to %.3f ms",
	     __FILE__, __LINE__, serv_info->exposure_time);

   /*
    * Set the default frame rate to be 50 Hz
    */
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) attempting to set the frame rate to %5.2f Hz",
	     __FILE__, __LINE__, DEFAULT_FRAME_RATE);
   if (setGuiderFrameRate(DEFAULT_FRAME_RATE) != PASS) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) unable to set frame rate to %5.2f Hz - exiting",
		__FILE__, __LINE__, DEFAULT_FRAME_RATE);
      exit(EXIT_FAILURE);
   }
   serv_info->frame_rate = DEFAULT_FRAME_RATE;
   
   cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
	     "(%s:%d) frame rate in camera set to %5.2f Hz",
	     __FILE__, __LINE__, DEFAULT_FRAME_RATE);

   /*
    * Create a linked list to hold client entries
    */
   serv_info->client_list = createList(cli_malloc);

   /* 
    * Set up the server parameters
    */
   if ((serv_info->raptor_serv = sockserv_create(RAPTOR_PORT)) == NULL) {
      cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		"(%s:%d) unable to establish a listening port for the"
		" server", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

   /* 
    * Cleanup camera and socket resources before exiting 
    */
   atexit(cleanup);

   /*
    * Register callback functions 
    */
   serv_info->raptor_serv->client_add_hook = clientAdd;
   serv_info->raptor_serv->client_del_hook = clientDelete;
   serv_info->raptor_serv->client_recv_hook = clientReceive;

   fprintf(stderr, "Ready to answer requests\n");

   /*
    * Handle termination and interrupt signals to close the shutter and exit
    */
   cli_signal(SIGTERM, cleanup);
   cli_signal(SIGINT, cleanup);

   /* 
    * Go through a loop processing any commands sent by the client.
    */
   for (;;) {

      cli_signal_block(SIGTERM);
      cli_signal_block(SIGINT);

      sockserv_run(serv_info->raptor_serv, SOCKSERV_POLL_INTERVAL);

      cli_signal_unblock(SIGTERM);
      cli_signal_unblock(SIGINT);
      
      /*
       * Determine if a request has been made to turn on video mode when it
       * was off
       */
      if ((last_video_on_state == FALSE) && (serv_info->video_on == TRUE)) {
      
	 /* 
	  * Try to open a handle to the device     
	  */
	 if (serv_info->pdv_p == NULL) {
	    if ((serv_info->pdv_p = pdv_open_channel(edt_devname, edt_unit, 
						     edt_channel)) == NULL) {
	       cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
			 "(%s:%d) pdv_open_channel request failed",
			 __FILE__, __LINE__);
	       serv_info->video_on = FALSE;
	       serv_info->pdv_p = NULL;
	       continue;
	    }
	 }
      
	 /*
	  * Get the width and height for the image
	  */
	 serv_info->image_width = pdv_get_width(serv_info->pdv_p);
	 serv_info->image_height = pdv_get_height(serv_info->pdv_p);

	 /*
	  * Make sure the height and width are valid
	  */
	 if (serv_info->image_width <= 1 && serv_info->image_height <= 1) {
	    cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		      "(%s:%d) image size incorrect (width=%d, height=%d)",
		      __FILE__, __LINE__, serv_info->image_width,
		      serv_info->image_height);
	    serv_info->video_on = FALSE;
	    pdv_close(serv_info->pdv_p);
	    serv_info->pdv_p = NULL;
	    continue;
	 }

	 /*
	  * The number of buffers is limited only by the amount of host 
	  * memory available, up to approximately 3.5 GBytes (or less, 
	  * depending on other OS use of the low 3.5 GB of memory.  Each 
	  * buffer has a certain amount of overhead, so setting a large
	  * number, event if the images are small, is not recommended.  Four
	  * is the recommended number; at any time, one buffer is being read
	  * in, one buffer is being read out, one is being set up for DMA, 
	  * and one is reserved in case of overlap.  Additional buffers may
	  * be necessary with very fast cameras.  32 will almost always 
	  * smooth out any problems with really fast cameras, and if the
	  * system can't keep up with 64 buffers allocated, there may be
	  * other problems.
	  */
	 if (pdv_multibuf(serv_info->pdv_p, 4) != 0) {
	    cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		      "(%s:%d) pdv_multibuf() call failed",
		      __FILE__, __LINE__);
	    serv_info->video_on = FALSE;
	    pdv_close(serv_info->pdv_p);
	    serv_info->pdv_p = NULL;
	    continue;
	 }

	 /*
	  * Set the timeout to block on pdv_wait_image
	  */
	 if (pdv_set_timeout(serv_info->pdv_p, 0) != 0) {
	    cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		      "(%s:%d) pdv_set_timeout() call failed",
		      __FILE__, __LINE__);
	    serv_info->video_on = FALSE;
	    pdv_close(serv_info->pdv_p);
	    serv_info->pdv_p = NULL;
	    continue;
	 }

	 /*
	  * Mark that the state indicates that video is on
	  */
	 last_video_on_state = TRUE;

#ifdef DEBUG
         /* 
          * Open a csv file where to write:
          * Index; time (ms); Xcmd (arcsec); Ycmd (arcsec);
          * Xmes (arcsec); Ymes (arcsec)
          */
         filePtr = fopen("/cfht/src/spirou/guider/raptorServ/Outputs.csv", "w");
         if (filePtr == NULL) {
            cfht_logv(CFHT_MAIN, CFHT_ERROR, "(%s:%d) Failed to open "
                      "output file", __FILE__, __LINE__);
             exit(EXIT_FAILURE);
         }
         /* Print out header to the csv file */
         if (fprintf(filePtr, "Index;Time(ms);Xcmd (arcsec);Ycmd (arcsec)"
                     ";Xmes (arcsec);Ymes (arcsec);DeltaX (arcsec);DeltaY (arcsec)\n") < 0) {
            cfht_logv(CFHT_MAIN, CFHT_ERROR, "(%s:%d) Failed to print to "
                      "output file", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
         }
#endif

      } /* End of video switch on */

      /*
       * Handle images for video mode
       */
      if (serv_info->video_on == TRUE) {

	 /*
	  * Start the acquisition of the next image
	  */
	 pdv_start_images(serv_info->pdv_p, 0);

	 /*
	  * Return the latest image
	  */
	 image_p = pdv_wait_image(serv_info->pdv_p);
	 
	 /*
	  *  Starting the centroid calculation
	  */
	 if (serv_info->guide_on == TRUE){
		 if (serv_info->fwhm_flag == 0){
		   calculatePointFWHM((unsigned short *)image_p, GUIDE_SIZE_X, GUIDE_SIZE_Y);
		   serv_info->fwhm_flag=1;

		   //fprintf(stderr, "%7.3f %7.3f \n",fwhm_x,fwhm_y);
		 }


	       /*
		*  Calculate the centroid
		*/
#ifndef SIM_STAR
	   //calculateCentroid((unsigned short *)image, GUIDE_SIZE_X, GUIDE_SIZE_Y,&xc, &yc);
	   calculateCentroidMPFIT((unsigned short *)image_p, GUIDE_SIZE_X, GUIDE_SIZE_Y,&xc, &yc);

	   serv_info->guide_xoff=xc;
	   serv_info->guide_yoff=yc;
#endif

	   if (serv_info->isu_on == TRUE){
               /* Checking ISU error status first of all */
               if (check_isu(&x_fault, &y_fault) == FAIL) {
	          cfht_logv(CFHT_MAIN, CFHT_WARN, "(%s:%d) Failed checking ISU",
		      __FILE__, __LINE__);
               }
               else {
                  if (x_fault) {
	             cfht_logv(CFHT_MAIN, CFHT_ERROR, "(%s:%d) Fatal error on "
                     "the fast guiding loop: ISU x axis is in error.  Relaunch "
                     "the fast guiding loop. An evolution of libisu to avoid "
                     "this fatal error should be considered", 
		      __FILE__, __LINE__);
                     if (xangle != 0) {
	                cfht_logv(CFHT_MAIN, CFHT_ERROR, "(%s:%d) \"Setup\" "
                        "position requested was %f arcsec in x",
                        __FILE__, __LINE__, xangle);
                     fflush(stdout);
                     }
                     exit(EXIT_FAILURE);
                  }
                  if (y_fault) {
	             cfht_logv(CFHT_MAIN, CFHT_ERROR, "(%s:%d) Fatal error on "
                     "the fast guiding loop: ISU y axis is in error.  Relaunch "
                     "the fast guiding loop. An evolution of libisu to avoid "
                     "this fatal error should be considered", 
		      __FILE__, __LINE__);
                     if (yangle != 0) {
	                cfht_logv(CFHT_MAIN, CFHT_ERROR, "(%s:%d) \"Setup\" "
                        "position requested was %f arcsec in y",
                        __FILE__, __LINE__, yangle);
                     fflush(stdout);
                     }
                     exit(EXIT_FAILURE);
                  }
               } /* End of ISU check */

               /* 
                * Retreiving current ISU position in mrad on the mechanism
                * ("true" position)
                */
               if (get_angles(&last_x_angle, &last_y_angle) != PASS) {
	          cfht_logv(CFHT_MAIN, CFHT_ERROR, "(%s:%d) Fatal error on "
                  "the fast guiding loop: failed getting ISU angles",
		      __FILE__, __LINE__);
                  exit(EXIT_FAILURE);
               }
#ifdef DEBUG
               /* 
                * At this place in the code next position is not yet 
                * defined so xangle is the last setup
                */
               gettimeofday(&this_time,&tz);

               /* time_spent is in ms */
               time_spent = ((double) this_time.tv_sec  * 1E3+
                      (double) this_time.tv_usec * 1E-3) -
                      ((double) last_time.tv_sec * 1E3  +
                      (double) last_time.tv_usec * 1E-3);
              /* 
               * Print out results to the csv file
               * Index; time (ms); Xcmd (arcsec); Ycmd (arcsec);
               * Xmes (arcsec); Ymes (arcsec);DeltaX (arcsec);DeltaY (arcsec)
               */
               index++;
               meas_setup_x_angle = last_x_angle;
               meas_setup_y_angle = last_y_angle;
               true_to_setup(&meas_setup_x_angle, &meas_setup_y_angle);
               if (fprintf(filePtr, "%d;%lf;%lf;%lf;%lf;%lf;%lf;%lf\n",
                   index, time_spent, (double)xangle, (double)yangle,
                   meas_setup_x_angle, meas_setup_y_angle, 
                   (double)xangle - meas_setup_x_angle,
                   (double)yangle - meas_setup_y_angle) < 0) {
	          cfht_logv(CFHT_MAIN, CFHT_ERROR, "(%s:%d) Failed to print to "
                  "output file", __FILE__, __LINE__);
                   exit(EXIT_FAILURE);
               }
               last_time.tv_sec = this_time.tv_sec;
               last_time.tv_usec = this_time.tv_usec;
#endif

#ifdef SIM_STAR
              /* 
               * Simulate a gaussian star distribution in x & y
               * sampleNormal function gives Gaussian distribution
               * with zero mean and a standard deviation of one
               * for ISU  tests we need ito simulate a seeing of 0,65 arcsec
               * so FWHM = 2.35486 sigma = 0.65 so sigma is 0.276 
               */
               sampleNormal(&xangle, &yangle);
               xangle = xangle * 0.276;
               yangle = yangle * 0.276;

	       serv_info->guide_xoff=0;
	       serv_info->guide_yoff=0;
#else
               /* Converting pixel values to angle in arcsec */
               xangle=(xc-(GUIDE_SIZE_X/2))* PIXSCALE;
               yangle=(yc-(GUIDE_SIZE_Y/2))* PIXSCALE;
#endif

               /* Filling in the thread_data structure */
               /* arg1 is FGL frequency in Hz. It is the frame_rate */
               if (serv_info->frame_rate != 0) {
                  thread_data.arg1 = (double)serv_info->frame_rate;
               }
               else {
                  thread_data.arg1 = (double)DEFAULT_FRAME_RATE;
               }
               /* arg2 is last x position in mrad */
               thread_data.arg2 = last_x_angle;
               /* arg3 is last y position in mrad */ 
               thread_data.arg3 = last_y_angle;
               next_x_angle = (double)xangle;
               next_y_angle = (double)yangle;
               setup_to_true(&next_x_angle, &next_y_angle);
               /* arg4 is next x position in mrad */
               thread_data.arg4 = next_x_angle;
               /* arg5 is next y position in mrad */
               thread_data.arg5 = next_y_angle;

              /*
	       *   Sending ISU corrections
	       */

               /* Create thread that will allow calling set_analog_slope */
               if (pthread_create(&set_analog_slope_t, NULL,
                          set_analog_slope_thread, (void *)&thread_data)) {
                  cfht_logv(CFHT_MAIN, CFHT_ERROR, "(%s:%d) Fatal error on "
                     "the fast guiding loop: failed creating generation thread",
                     __FILE__, __LINE__);
                  fflush(stdout);
                  exit(EXIT_FAILURE);
               }
               if (pthread_detach(set_analog_slope_t) != PASS) {
                  cfht_logv(CFHT_MAIN, CFHT_WARN,
                         "(%s:%d) %s: Failed detaching separate thread !",
                        __FILE__, __LINE__, __FUNCTION__);
                  fflush(stdout);
                  exit(EXIT_FAILURE);
               }
              /* OR... send corrections without thread */
/*
              if (setup_isu(next_x_angle, next_y_angle) != PASS) {
                  cfht_logv(CFHT_MAIN, CFHT_WARN,
                         "(%s:%d) %s: Error setting isu corrections",
                        __FILE__, __LINE__, __FUNCTION__);
              }
*/
	   } /* End of isu correction loop */
	 } /* End of guide loop */

	 /*
	  * Check if any timeouts occured since the last image
	  */
	 timeouts = pdv_timeouts(serv_info->pdv_p);
	 if (timeouts > last_timeouts) {
	    cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		      "(%s:%d) received %d timeouts since the last image",
		      __FILE__, __LINE__, timeouts - last_timeouts);
	 }
	 
	 /*
	  * Create a FITS image from the pixel data and send it to stdout
	  */
	 if (writeFITSImage(image_p) != PASS) {
	    cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		      "(%s:%d) unable to create FITS file and write it to"
		      " STDOUT", __FILE__, __LINE__);
	 }
      }

      /*
       * If the video mode is turned off, change the last state to reflect
       * this.
       */
      if (serv_info->video_on == FALSE) {
         if (last_video_on_state == TRUE) {
#ifdef DEBUG
            /* Close the debug csv file */
            if (fclose(filePtr) != 0) {
               cfht_logv(CFHT_MAIN, CFHT_WARN, "(%s:%d) Failed to close "
                  "output file", __FILE__, __LINE__);
            }
#endif
	    last_video_on_state = FALSE;
         }
      }
   }

   exit(EXIT_SUCCESS);
}
