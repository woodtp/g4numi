#!/bin/bash

# use the bootstrap script to build a source code distribution

# bootstrap.sh <product_directory> 

usage()
{
   echo "USAGE: `basename ${0}` <product_dir>"
   echo "       `basename ${0}` installs the source code distribution and creates the source code tarball"
}


get_this_dir() 
{
    ( cd / ; /bin/pwd -P ) >/dev/null 2>&1
    if (( $? == 0 )); then
      pwd_P_arg="-P"
    fi
    reldir=`dirname ${0}`
    thisdir=`cd ${reldir} && /bin/pwd ${pwd_P_arg}`
}

productdir=${1}

if [ -z ${productdir} ]
then
   echo "ERROR: please specify the product directory"
   usage
   exit 1
fi

package=novasoft
pkgver=devel
pkgdotver=`echo ${pkgver} | sed -e 's/_/./g' | sed -e 's/^v//'`
sourceurl=""

tarballname=${package}-${pkgdotver}-source.tar.bz2
pkgdir=${productdir}/${package}/${pkgver}

get_this_dir

# remove any existing distribution
if [ -d ${pkgdir} ]
then
    echo "removing previous install of ${package} ${pkgver}"
    rm -rf ${pkgdir}
    rm -rf ${pkgdir}.version
fi
# now build the new package source code distribution
mkdir -p ${pkgdir}
if [ ! -d ${pkgdir} ]
then
   echo "ERROR: unable to create ${pkgdir}"
   exit 1
fi
# install the build script
if [ ! -e ${thisdir}/build_${package}.sh ]
then
   echo "ERROR: cannot find ${thisdir}/build_${package}.sh"
   exit 1
else
   cp -p ${thisdir}/build_${package}.sh ${pkgdir}/
fi
# install the autobuild script
if [ ! -e ${thisdir}/autobuild.sh ]
then
   echo "ERROR: cannot find ${thisdir}/autobuild.sh"
   exit 1
else
   cp -p ${thisdir}/autobuild.sh ${pkgdir}/
fi

set -x
# get the source code 
cd ${pkgdir}

svn export http://cdcvs.fnal.gov/subversion/novaart.pkgs.svn/trunk cet

# make the source code tarball
cd ${productdir}
tar cjf ${tarballname} \
   ${package}/${pkgver}/*.sh \
   ${package}/${pkgver}/cet
cd ${thisdir}

set +x


exit 0

