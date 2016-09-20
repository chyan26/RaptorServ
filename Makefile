# `Makefile' - Use with GNU make to build ssGetenv
#
#   This file is part of version 1 of the environment server shell client.
#   Read the `License' file for terms of use and distribution.
#   Copyright 2002, Canada-France-Hawaii Telescope, daprog@cfht.hawaii.edu.
#
# ___This header automatically generated from `Index'. Do not edit it here!___

include ../Make.Common

#CCWARN += $(WERROR)
LOCALLIBS += libisu.a
LOCALLIBS += /cfht/src/spirou/guider/powerdaq-3.6.24/lib/libpowerdaq32.so.1.0
$(EXECNAME) $(EXECNAME)-pure: $(OBJS)   
EXTRA_CCLINK += -lfh -lcli -lcfht -lm -lsgc -lpthread -lpdv -ldl -lmpfit -lcfitsio -lisu -lpowerdaq32 -lsockio -lssapi -lss
CCINCS += -I/cfht/include/isu/
include ../Make.Common

# Dependencies by Make.Common $Revision: 2.17 $
