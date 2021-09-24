#############################################################################
# Makefile for building: simul_ecg
# Generated by qmake (3.1) (Qt 6.0.0)
# Project:  simul_ecg.pro
# Template: app
# Command: /home/neil/Qt/6.0.0/gcc_64/bin/qmake -o Makefile simul_ecg.pro -spec linux-g++ CONFIG+=qtquickcompiler
#############################################################################

MAKEFILE      = Makefile

EQ            = =

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = 
CFLAGS        = -pipe -O2 -Wall -Wextra -fPIC $(DEFINES)
CXXFLAGS      = -pipe -O3 -DNDEBUG -O2 -std=gnu++1z -Wall -Wextra -fPIC $(DEFINES)
INCPATH       = -I. -I../../../Qt/6.0.0/gcc_64/mkspecs/linux-g++
QMAKE         = /home/neil/Qt/6.0.0/gcc_64/bin/qmake
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
QINSTALL      = /home/neil/Qt/6.0.0/gcc_64/bin/qmake -install qinstall
QINSTALL_PROGRAM = /home/neil/Qt/6.0.0/gcc_64/bin/qmake -install qinstall -exe
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
DISTNAME      = simul_ecg1.0.0
DISTDIR = /home/neil/Workspace/CPP/simul_ecg/.tmp/simul_ecg1.0.0
LINK          = g++
LFLAGS        = -Wl,-O1
LIBS          = $(SUBLIBS) -L/usr/local/lib -lstdc++fs -lpthread   
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = strip

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = main.cc \
		opt.c 
OBJECTS       = main.o \
		opt.o
DIST          = 00README \
		index.html \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/spec_pre.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/unix.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/linux.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/sanitize.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/gcc-base.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/gcc-base-unix.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/g++-base.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/g++-unix.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/qconfig.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_ext_libpng.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_concurrent.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_concurrent_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_core.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_core_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_core_qobject_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_dbus.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_dbus_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_designer.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_designer_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_designercomponents_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_devicediscovery_support_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_eglfs_kms_support_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_eglfsdeviceintegration_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_fb_support_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_gui.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_gui_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_help.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_help_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_input_support_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_kms_support_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_linguist.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_linguist_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_network.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_network_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_opengl.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_opengl_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_openglwidgets.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_openglwidgets_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_packetprotocol_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_printsupport.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_printsupport_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qml.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qml_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlcompiler_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmldebug_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmldevtools_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlmodels.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlmodels_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmltest.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmltest_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlworkerscript.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlworkerscript_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quick.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quick_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2impl.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2impl_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickparticles_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickshapes_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quicktemplates2.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quicktemplates2_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickwidgets.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickwidgets_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_sql.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_sql_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svg.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svg_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svgwidgets.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svgwidgets_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_testlib.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_testlib_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_tools_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_uiplugin.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_uitools.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_uitools_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_waylandclient.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_waylandclient_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_widgets.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_widgets_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_xml.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_xml_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/qt_functions.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/qt_config.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/linux-g++/qmake.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/spec_post.prf \
		.qmake.stash \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/exclusive_builds.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/toolchain.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/default_pre.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/resolve_config.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/default_post.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/qtquickcompiler.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/warn_on.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/qmake_use.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/file_copies.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/testcase_targets.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/exceptions.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/yacc.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/lex.prf \
		simul_ecg.pro dfour1.h \
		dir_utils.hh \
		opt.h \
		options.hh \
		rand1.h \
		single_include/csv.hpp main.cc \
		opt.c
QMAKE_TARGET  = simul_ecg
DESTDIR       = 
TARGET        = simul_ecg


first: all
####### Build rules

simul_ecg:  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: simul_ecg.pro ../../../Qt/6.0.0/gcc_64/mkspecs/linux-g++/qmake.conf ../../../Qt/6.0.0/gcc_64/mkspecs/features/spec_pre.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/unix.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/linux.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/sanitize.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/gcc-base.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/gcc-base-unix.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/g++-base.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/common/g++-unix.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/qconfig.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_ext_libpng.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_concurrent.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_concurrent_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_core.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_core_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_core_qobject_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_dbus.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_dbus_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_designer.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_designer_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_designercomponents_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_devicediscovery_support_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_eglfs_kms_support_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_eglfsdeviceintegration_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_fb_support_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_gui.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_gui_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_help.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_help_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_input_support_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_kms_support_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_linguist.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_linguist_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_network.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_network_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_opengl.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_opengl_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_openglwidgets.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_openglwidgets_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_packetprotocol_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_printsupport.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_printsupport_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qml.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qml_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlcompiler_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmldebug_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmldevtools_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlmodels.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlmodels_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmltest.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmltest_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlworkerscript.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlworkerscript_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quick.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quick_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2impl.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2impl_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickparticles_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickshapes_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quicktemplates2.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quicktemplates2_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickwidgets.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickwidgets_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_sql.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_sql_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svg.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svg_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svgwidgets.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svgwidgets_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_testlib.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_testlib_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_tools_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_uiplugin.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_uitools.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_uitools_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_waylandclient.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_waylandclient_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_widgets.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_widgets_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_xml.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_xml_private.pri \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/qt_functions.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/qt_config.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/linux-g++/qmake.conf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/spec_post.prf \
		.qmake.stash \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/exclusive_builds.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/toolchain.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/default_pre.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/resolve_config.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/default_post.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/qtquickcompiler.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/warn_on.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/qmake_use.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/file_copies.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/testcase_targets.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/exceptions.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/yacc.prf \
		../../../Qt/6.0.0/gcc_64/mkspecs/features/lex.prf \
		simul_ecg.pro
	$(QMAKE) -o Makefile simul_ecg.pro -spec linux-g++ CONFIG+=qtquickcompiler
../../../Qt/6.0.0/gcc_64/mkspecs/features/spec_pre.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/common/unix.conf:
../../../Qt/6.0.0/gcc_64/mkspecs/common/linux.conf:
../../../Qt/6.0.0/gcc_64/mkspecs/common/sanitize.conf:
../../../Qt/6.0.0/gcc_64/mkspecs/common/gcc-base.conf:
../../../Qt/6.0.0/gcc_64/mkspecs/common/gcc-base-unix.conf:
../../../Qt/6.0.0/gcc_64/mkspecs/common/g++-base.conf:
../../../Qt/6.0.0/gcc_64/mkspecs/common/g++-unix.conf:
../../../Qt/6.0.0/gcc_64/mkspecs/qconfig.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_ext_libpng.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_concurrent.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_concurrent_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_core.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_core_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_core_qobject_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_dbus.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_dbus_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_designer.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_designer_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_designercomponents_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_devicediscovery_support_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_eglfs_kms_support_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_eglfsdeviceintegration_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_fb_support_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_gui.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_gui_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_help.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_help_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_input_support_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_kms_support_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_linguist.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_linguist_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_network.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_network_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_opengl.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_opengl_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_openglwidgets.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_openglwidgets_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_packetprotocol_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_printsupport.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_printsupport_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qml.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qml_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlcompiler_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmldebug_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmldevtools_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlmodels.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlmodels_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmltest.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmltest_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlworkerscript.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_qmlworkerscript_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quick.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quick_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2impl.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickcontrols2impl_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickparticles_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickshapes_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quicktemplates2.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quicktemplates2_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickwidgets.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_quickwidgets_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_sql.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_sql_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svg.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svg_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svgwidgets.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_svgwidgets_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_testlib.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_testlib_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_tools_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_uiplugin.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_uitools.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_uitools_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_waylandclient.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_waylandclient_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_widgets.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_widgets_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_xml.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/modules/qt_lib_xml_private.pri:
../../../Qt/6.0.0/gcc_64/mkspecs/features/qt_functions.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/qt_config.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/linux-g++/qmake.conf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/spec_post.prf:
.qmake.stash:
../../../Qt/6.0.0/gcc_64/mkspecs/features/exclusive_builds.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/toolchain.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/default_pre.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/resolve_config.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/default_post.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/qtquickcompiler.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/warn_on.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/qmake_use.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/file_copies.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/testcase_targets.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/exceptions.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/yacc.prf:
../../../Qt/6.0.0/gcc_64/mkspecs/features/lex.prf:
simul_ecg.pro:
qmake: FORCE
	@$(QMAKE) -o Makefile simul_ecg.pro -spec linux-g++ CONFIG+=qtquickcompiler

qmake_all: FORCE


all: Makefile simul_ecg

dist: distdir FORCE
	(cd `dirname $(DISTDIR)` && $(TAR) $(DISTNAME).tar $(DISTNAME) && $(COMPRESS) $(DISTNAME).tar) && $(MOVE) `dirname $(DISTDIR)`/$(DISTNAME).tar.gz . && $(DEL_FILE) -r $(DISTDIR)

distdir: FORCE
	@test -d $(DISTDIR) || mkdir -p $(DISTDIR)
	$(COPY_FILE) --parents $(DIST) $(DISTDIR)/


clean: compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


distclean: clean 
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) .qmake.stash
	-$(DEL_FILE) Makefile


####### Sub-libraries

check: first

benchmark: first

compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: 

####### Compile

main.o: main.cc dir_utils.hh \
		single_include/csv.hpp \
		opt.h \
		options.hh \
		rand1.h \
		dfour1.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cc

opt.o: opt.c 
	$(CC) -c $(CFLAGS) $(INCPATH) -o opt.o opt.c

####### Install

install:  FORCE

uninstall:  FORCE

FORCE:
