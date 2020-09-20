#!/bin/bash

# WEIS CCD Tool Installation Script
# Written by Yong Hoon Lee (ylee196@illinois.edu)
# This script is for Ubuntu Linux and tested with 20.04 LTS

# Prerequisites:
# - Anaconda 3 with Python 3.8 or higher (environment will be created)
# - Intel Parallel Studio XE (Fortran, C, C++)
# - Intel Math Kernel Library

# If root privilege is needed (for now, it is not required)
#if [ $(id -u) -ne 0 ]; then
#	echo 'System root privilege is needed to run this script'
#	exit 1
#fi

# Parameters
inteldir='/opt/intel'

# Checking Anaconda

echo '=============================='
echo 'Checking anaconda installation'
echo '------------------------------'
condainitscript=$(cat ~/.bashrc | sed -n '/__conda_setup=/,/unset __conda_setup/p')
if [ -z "$condainitscript" ]; then
	echo 'ERROR: Anaconda should be installed and initialized.'
	exit 1
else
	condainstalldirectory=$(which conda)
	if [ -z "$condainstalldirectory" ]; then
		echo 'ERROR: Anaconda path is not recongnized.'
		exit 1
	else
		echo 'Anaconda path is recognized as: '$condainstalldirectory'.'
		if [[ $condainstalldirectory == $HOME* ]]; then
			echo 'Anaconda is installed in your home directory.'
		else
			echo 'ERROR: Anaconda should be installed in your home directory.'
			exit 1
		fi
	fi
fi

# Checking Intel Parallel Studio XE (works with 2020)

echo '=============================================='
echo 'Checking Intel Parallel Studio XE installation'
echo '----------------------------------------------'
if [ -d $inteldir ]; then
	echo 'Intel software packages are located in '$inteldir'.'
else
	echo 'Unable to locate Intel software packages at '$inteldir'.'
	exit 1
fi

intelcompilerscript=$(cat ~/.bashrc | grep 'psxevars.sh')
if [ -z "$intelcompilerscript" ]; then
	echo 'Unable to find environment configuration of Intel Parallel Studio XE in your .bashrc file.'
	intelpsxelnk=$(find $inteldir -maxdepth 1 -type l -name 'parallel_studio_xe*')
	intelpsxesource=${intelpsxelnk}/bin/psxevars.sh
	if [ -f "$intelpsxesource" ]; then
		echo 'psxevars.sh found.'
		if [ -f ~/.bashrc.weis_backup_psxesource ]; then
			rm ~/.bashrc.weis_backup_psxesource
		fi
		cp ~/.bashrc ~/.bashrc.weis_backup_psxesource
		sed -i '/# >>> conda initialize >>>/i # Intel Parallel Studio XE' ~/.bashrc
		sed -i '/# >>> conda initialize >>>/i source '$intelpsxesource ~/.bashrc
		sed -i '/# >>> conda initialize >>>/i \\' ~/.bashrc
		echo 'Intel Parallel Studio XE environment is configured in your .bashrc file.'
		echo 'ATTENTION: Please restart bash terminal and run this script again.'
		exit 1
	else
		echo 'ERROR: Unable to find psxevars.sh. Please install Intel Parallel Studion XE.'
		exit 1
	fi
else
	echo 'Found environment configuration of Intel Parallel Studio XE from your .bashrc file.'
fi

# Checking Intel Compiler Suite

echo '=========================================='
echo 'Checking Intel Compiler Suite Installation'
echo '------------------------------------------'
flg=0
icc=$(which icc)
if [ -z "$icc" ]; then
	echo 'ERROR: Intel C compiler (icc) not found.'
	exit 1
else
	echo 'Intel C compiler (icc) is located at '$icc'.'
	if [ -z "$CC" ]; then
		echo 'Unable to find CC parameter.'
		if [ -f ~/.bashrc.weis_backup_ccparam ]; then
			rm ~/.bashrc.weis_backup_ccparam
		fi
		cp ~/.bashrc ~/.bashrc.weis_backup_ccparam
		sed -i '/# >>> conda initialize >>>/i # Compiler Parameters' ~/.bashrc
		sed -i '/# >>> conda initialize >>>/i export CC='$icc ~/.bashrc
		echo 'CC parameter is added for Intel C compiler.'
		flg=1
	else
		echo 'Checking environment variable CC'
		if [[ "$CC" == *icc ]]; then
			echo 'CC parameter is set to Intel C compiler.'
		elif [[ "$CC" == *gcc ]]; then
			echo 'WARNING: CC parameter is set to GNU C compiler.'
		else
			echo 'WARNING: CC parameter is set to unknown C compiler.'
		fi
	fi
fi
icpc=$(which icpc)
if [ -z "$icpc" ]; then
	echo 'ERROR: Intel C++ compiler (icpc) not found.'
	exit 1
else
	echo 'Intel C++ compiler (icpc) is located at '$icpc'.'
	if [ -z "$CXX" ]; then
		echo 'Unable to find CXX parameter.'
		if [ -f ~/.bashrc.weis_backup_cxxparam ]; then
			rm ~/.bashrc.weis_backup_cxxparam
		fi
		cp ~/.bashrc ~/.bashrc.weis_backup_cxxparam
		if [[ "$flg" == 0 ]]; then
			sed -i '/# >>> conda initialize >>>/i # Compiler Parameters' ~/.bashrc
		fi
		sed -i '/# >>> conda initialize >>>/i export CXX='$icpc ~/.bashrc
		echo 'CXX parameter is added for Intel C++ compiler.'
		flg=1
	else
		echo 'Checking environment variable CXX'
		if [[ "$CXX" == *icpc ]]; then
			echo 'CXX parameter is set to Intel C++ compiler.'
		elif [[ "$CXX" == *g++ ]]; then
			echo 'WARNING: CXX parameter is set to GNU C++ compiler.'
		else
			echo 'WARNING: CXX parameter is set to unknown C++ compiler.'
		fi
	fi
fi
ifort=$(which ifort)
if [ -z "$ifort" ]; then
	echo 'ERROR: Intel Fortran compiler (ifort) not found.'
	exit 1
else
	echo 'Intel Fortran compiler (ifort) is located at '$ifort'.'
	if [ -z "$FC" ]; then
		echo 'Unable to find FC parameter.'
		if [ -f ~/.bashrc.weis_backup_fcparam ]; then
			rm ~/.bashrc.weis_backup_fcparam
		fi
		cp ~/.bashrc ~/.bashrc.weis_backup_fcparam
		if [[ "$flg" == 0 ]]; then
			sed -i '/# >>> conda initialize >>>/i # Compiler Parameters' ~/.bashrc
		fi
		sed -i '/# >>> conda initialize >>>/i export FC='$ifort ~/.bashrc
		echo 'FC parameter is added for Intel Fortran compiler.'
		flg=1
	else
		echo 'Checking environment variable FC'
		if [[ "$FC" == *ifort ]]; then
			echo 'FC parameter is set to Intel Fortran compiler.'
		elif [[ "$FC" == *gfortran ]]; then
			echo 'WARNING: FC parameter is set to GNU Fortran compiler.'
		else
			echo 'WARNING: FC parameter is set to unknown Fortran compiler.'
		fi
	fi
fi
if [[ "$flg" == 1 ]]; then
	sed -i '/# >>> conda initialize >>>/i \\' ~/.bashrc
	echo 'ATTENTION: Please restart bash terminal and run this script again.'
	exit 1
fi

# Checking Intel Math Kernel Library

echo '==============================================='
echo 'Checking Intel Math Kernel Library Installation'
echo '-----------------------------------------------'
if [ -z "$MKLROOT" ]; then
	echo 'ERROR: Intel Math Kernel Library (MKL) not found.'
	exit 1
else
	echo 'Intel Math Kernel Library (MKL) is located at '$MKLROOT'.'
	if [ -d $__mkl_lib_dir ]; then
		echo 'Intel Math Kernel Library (MKL) library path is recognized.'
	else
		echo 'ERROR: Missing required internal directories. Please install Intel Math Kernel Library (MKL)'
		exit 1
	fi
fi


