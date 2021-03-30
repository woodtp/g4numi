#!/bin/bash
#

# build_nova-daq-offline.sh <product_dir> <e4|e5> <debug|opt|prof> [tar] [source]

usage()
{
   echo "USAGE: `basename ${0}` <product_dir> <e2|e4> <debug|opt|prof> [tar]"
}

# -------------------------------------------------------------------
# shared boilerplate
# -------------------------------------------------------------------

get_this_dir() 
{
    ( cd / ; /bin/pwd -P ) >/dev/null 2>&1
    if (( $? == 0 )); then
      pwd_P_arg="-P"
    fi
    reldir=`dirname ${0}`
    thisdir=`cd ${reldir} && /bin/pwd ${pwd_P_arg}`
}

make_tarball()
{
    echo ""
    set -x
    cd ${product_dir}
    tarname=${package}-${pkgdotver}-${OS}-${plat}-${qualdir}.tar.bz2
    tar cjf ${tarname} \
       ${package}/${pkgver}.version/${flvr}_${qualver} \
       ${package}/${pkgver}/ups \
       ${package}/${pkgver}/source \
       ${package}/${pkgver}/include \
       ${package}/${pkgver}/job \
       ${package}/${pkgver}/pwd \
       ${package}/${pkgver}/xml \
       ${package}/${pkgver}/data \
       ${package}/${pkgver}/${OS}.${plat}.${dqual}.${extraqual}
    cd ${thisdir}
    set +x
    exit 0
}

make_source_tarball()
{
    set -x
    cd ${product_dir}
    tar cjf ${package}-${pkgdotver}-source.tar.bz2 \
       ${package}/${pkgver}/*.sh \
       ${package}/${pkgver}/cet
    cd ${thisdir}
    set +x
}

set_qualifiers()
{
 
    if [ -z ${basequal} ]
    then
       echo "ERROR: please specify the qualifiers"
       usage
       exit 1
    fi

    fullqual=${basequal}:${extraqual}
    qualdir=`echo ${basequal} | sed -e 's/:/-/g'`-${extraqual}
    dqual=`echo ${basequal} | sed -e s%:%.%g`
    if [ "${extraqual}" = "opt" ]
    then
	qualver=`echo ${basequal} | sed -e 's/:/_/g'`_${extraqual}
    elif [ "${extraqual}" = "debug" ]
    then
	qualver=`echo ${basequal} | sed -e 's/:/_/g'`_${extraqual}
    elif [ "${extraqual}" = "prof" ]
    then
	qualver=`echo ${basequal} | sed -e 's/:/_/g'`_${extraqual}
    else
       echo "ERROR: please specify debug, opt, or prof"
       usage
       exit 1
    fi
}

check_with_lsb_release ()
{
   plat=`uname -p`
   if [ `lsb_release -d | cut -f2 | cut  -f1 -d" "` = "Ubuntu" ]
   then
      OS1=ub
      OSnum=`lsb_release -r | cut -f2 | cut -f1 -d"."`
      OS=${OS1}${OSnum}
      plat=`uname -m`
   elif [ ${OS1} = "Linux" ]
   then
      OSnum=`lsb_release -r | cut -f2 | cut -f1 -d"."`
      OS=${OS1}${OSnum}
      # Scientific Linux - slf should work
      if [ `lsb_release -d | cut  -f3 -d" "` = "SL" ]
      then
	 OS=slf${OSnum}
      # Scientific Linux Fermi
      elif [ `lsb_release -d | cut  -f3 -d" "` = "SLF" ]
      then
	 OS=slf${OSnum}
      #
      elif [ `lsb_release -i | cut -f2` = "ScientificFermi" ]
      then
	 OS=slf${OSnum}
      # pretend that SL6 is the same as SLF6
      elif [ `lsb_release -i | cut -f2` = "Scientific" ]
      then
	 OS=slf${OSnum}
      # pretend that CentOS is SLF
      elif [ `lsb_release -i | cut -f2` = "CentOS" ]
      then
	 OS=slf${OSnum}
      # pretend that RedHatEnterpriseServer is SLF
      elif [ `lsb_release -i | cut -f2` = "RedHatEnterpriseServer" ]
      then
	 OS=slf${OSnum}
      # Scientific Linux CERN
      elif [ `lsb_release -d | cut  -f4 -d" "` = "SLC" ]
      then
	 OS=slc${OSnum}
      elif [ `lsb_release -d | cut  -f4 -d" "` = "LTS" ]
      then
	 OS=slf${OSnum}
      # unrecognized - pretend that this is SLF
      else
	 OS=slf${OSnum}
      fi
   fi
}

default_names ()
{
      OSnum1=`uname -r | cut -f1 -d"."`
      OSnum2=`uname -r | cut -f2 -d"."`
      OS=${OS1}${OSnum1}${OSnum2}
      plat=`uname -p`
}

check_linux ()
{
   plat=`uname -p`
   if [ -e /etc/system-release-cpe ]
   then
      OSnum=`cat /etc/system-release-cpe | cut -f5 -d":" |  cut -f1 -d"."`
      if [ `cat /etc/redhat-release | cut  -f4 -d" "` = "SLC" ]
      then
         OS=slc${OSnum}
      else
         OS=slf${OSnum}
      fi
   elif [ -e /etc/redhat-release ]
   then
      if [ `cat /etc/redhat-release | cut  -f4 -d" "` = "SLC" ]
      then
         OS=slc${OSnum}
         OSnum=`cat /etc/redhat-release | cut -f6 -d" " |  cut -f1 -d"."`
      elif [ `cat /etc/redhat-release | cut  -f3 -d" "` = "Fermi" ]
      then
         OS=slf${OSnum}
         OSnum=`cat /etc/redhat-release | cut -f5 -d" " |  cut -f1 -d"."`
      elif [ `cat /etc/redhat-release | cut  -f3 -d" "` = "SLF" ]
      then
         OS=slf${OSnum}
         OSnum=`cat /etc/redhat-release | cut -f5 -d" " |  cut -f1 -d"."`
      else
         OS=slf${OSnum}
         OSnum=`cat /etc/redhat-release | cut -f5 -d" " |  cut -f1 -d"."`
      fi
   else
      default_names
   fi
}

get_platform_name()
{
    # make sure we can use the setup alias
    if [ -z ${UPS_DIR} ]
    then
       echo "ERROR: please setup ups"
       exit 1
    fi
    source `${UPS_DIR}/bin/ups setup ${SETUP_UPS}`

    OS1=`uname`
    flvr=`ups flavor`
    if [ ${OS1} = "Darwin" ]
    then
	flvr=`ups flavor -2`
	plat=`uname -m`
	OSnum=`uname -r | cut -f1 -d"."`
	#OS=${OS1}${OSnum}
	#OS=osx${macos1}${macos2}
        OS=d${OSnum}
    elif [ "${flvr0}" = "Linuxppc" ]
    then
       OS=slf5
       plat=`uname -p`
    elif [ -e /usr/bin/lsb_release ]
    then
       check_with_lsb_release
    elif [ ${OS1} = "Linux" ]
    then
       check_linux
    else
       default_names
    fi

    if [ -z ${flvr} ]
    then
       echo "ERROR: could not determine flavor for this platform"
       echo "       something is very wrong"
       exit 1
    fi

    pkgdir=${product_dir}/${package}/${pkgver}/${OS}.${plat}.${dqual}.${extraqual}
    sourcedir=${product_dir}/${package}/${pkgver}/cet
    builddir=${product_dir}/${package}/${pkgver}/build/${flvr}.${dqual}.${extraqual}
}

select_gcc()
{
    if [ "${basequal}" = "e2" ] ||  [ "${basequal}" = "gcc47" ]
    then
       gccdotver=4.7
       gccver=v4_7_1
    elif [ "${basequal}" = "e4" ] ||  [ "${basequal}" = "gcc48" ]
    then
       gccdotver=4.8
       gccver=v4_8_1
    elif [ "${basequal}" = "gcc46" ]
    then
       gccdotver=4.6
       gccver=v4_6_1
    fi
}

check_gcc()
{
    gccshortver=`gcc --version | head -1 | cut -f3 -d" " | cut -f1-2 -d"."`
    gccfullver=`gcc --version | head -1 | cut -f3 -d" "`
    if [ ${gccshortver} = ${gccdotver} ]
    then
      echo " OK - gcc version is ${gccfullver}"
    else
      echo "ERROR: wrong gcc version ${gccfullver}"
      exit 1
    fi
}

setup_gcc()
{
    setup gcc ${gccver}
    check_gcc
}

declare_product()
{
    echo "Declaring ${package} ${pkgver} for ${flvr} and ${fullqual}"
    if [ ${OS1} = "Darwin" ]
    then
       ups declare ${package} ${pkgver} -r ${package}/${pkgver} -2 -m ${package}.table -q ${fullqual} -z ${product_dir}
    else 
       ups declare ${package} ${pkgver} -r ${package}/${pkgver} -4 -m ${package}.table -q ${fullqual} -z ${product_dir}
    fi
}

# -------------------------------------------------------------------
# start processing
# -------------------------------------------------------------------

product_dir=${1}
basequal=${2}
extraqual=${3}
maketar=${4}
sourcetar=${5}

# -------------------------------------------------------------------
# package name and version
# -------------------------------------------------------------------

package=novasoft
pkgver=devel
pkgdotver=`echo ${pkgver} | sed -e 's/_/./g' | sed -e 's/^v//'`

get_this_dir

if [ "${sourcetar}" = "source" ]
then
  make_source_tarball
fi

set_qualifiers

get_platform_name

select_gcc

if [ "${maketar}" = "tar" ] && [ -d ${pkgdir}/lib ]
then
   make_tarball
fi

echo "building ${package} for ${OS}-${plat}-${qualdir} (flavor ${flvr})"

# -------------------------------------------------------------------
# package dependent stuff goes here
# -------------------------------------------------------------------
if [ "${extraqual}" = "opt" ]
then
    extra_opt="-o ${basequal}"
elif [ "${extraqual}" = "debug" ]
then
    extra_opt="-d ${basequal}"
elif [ "${extraqual}" = "prof" ]
then
    extra_opt="-p ${basequal}"
fi

# start with clean build subdirectories
if [ -d ${builddir} ]
then
   rm -rf ${builddir}
fi

echo "creating  ${builddir}"
mkdir -p ${builddir}
if [ ! -d ${builddir} ]
then
   echo "ERROR: failed to create ${builddir}"
   echo "       Something is very wrong"
   exit 1
fi

cd ${builddir}
source ${sourcedir}/ups/setup_for_development ${extra_opt}
set -x

env CC=gcc CXX=g++ FC=gfortran cmake -DCMAKE_INSTALL_PREFIX=${product_dir} -DCMAKE_BUILD_TYPE=$CETPKG_TYPE $CETPKG_SOURCE
make -j4
make ARGS=-j4 test
make install

set +x

if [ ! -d ${pkgdir}/lib ]
then
   echo "ERROR: failed to create ${pkgdir}/lib"
   echo "       Something is very wrong"
   exit 1
fi

# -------------------------------------------------------------------
# common bottom stuff
# -------------------------------------------------------------------

# this should not complain
echo "Finished building ${package} ${pkgver}"
setup ${package} ${pkgver} -q ${fullqual} -z ${product_dir}
echo "${package} is installed at ${NOVASOFT_FQ_DIR}"

# this must be last because it will exit
if [ "${maketar}" = "tar" ] && [ -d ${pkgdir}/lib ]
then
   make_tarball
fi

exit 0
