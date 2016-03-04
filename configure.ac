#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([bs_sim], [1.0], [eloi.casals@cnag.crg.eu])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([1.12 foreign no-define])
AC_CONFIG_SRCDIR([src/bs_sim.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [powf])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main])
AC_CHECK_LIB([pthread], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_CHECK_FUNCS([strchr])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_CONFIG_SUBDIRS([htslib-1.3])
AC_OUTPUT
