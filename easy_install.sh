#!/bin/bash -Ee

prompt_download() {
	echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	echo "About to install $1"
	echo "I will download from $2"
	read -p "OK? (Yes/abort)" yn
	case $yn in
		[Aa]* ) exit;;	
	esac
}

confirmed_download() {
	prompt_download $1 $2
	curl -L ${2} > tmp.tgz
	tar -xvzpf tmp.tgz
	rm -f tmp.tgz
}

install_zlib () {
	if [ -x "$(which apt-get)" ]; then
		apt-get install libgcrypt11-dev zlib1g-dev
	elif [ -x "$(which brew)" ]; then
		brew install zlib
	else 
		rm -Rf zlib*
		confirmed_download 'Zlib' 'http://zlib.net/zlib-1.2.8.tar.gz'
		cd zlib*
			./configure
			sudo make install
		cd ..
	fi
}

install_curl() {
	if [ -x "$(which apt-get)" ]; then
		apt-get install curl
	elif [ -x "$(which brew)" ]; then
		brew install curl
	else
		rm -Rf curl*
		confirmed_download 'curl' 'https://curl.haxx.se/download/curl-7.51.0.tar.gz'
		echo "Installing curl..."
		cd curl*
			./configure
			make
			sudo make install
		cd ..
	fi
}

install_gsl() {
	if [ -x "$(which apt-get)" ]; then
		apt-get install libgsl0ldbl
	elif [ -x "$(which brew)" ]; then
		brew install gsl
	else
		rm -Rf gsl*
		confirmed_download 'GSL' 'ftp://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz'
		echo "Installing GNU Scientific library..."
		cd gsl*
			./configure
			make
			sudo make install
		cd ..
	fi
}

install_brew() {
	prompt_download "HomeBrew" "https://raw.github.com/Homebrew/homebrew/go/install"
	ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/homebrew/go/install)"
}

install_git() {
	if [ -x "$(which apt-get)" ]; then
		apt-get install git
	elif [ -x "$(which brew)" ]; then
		brew install git
	else
		prompt_download "Git" "https://github.com/git/git/archive/master.zip"
		curl -L "https://github.com/git/git/archive/master.zip" > tmp.zip
		unzip -j master.zip -d git
		rm master.zip
		cd git*
			make prefix=/usr/local all
			make prefix=/usr/local install
		cd ..
	fi
	export PATH=$PATH
}

install_htslib() {
	prompt_download 'HTSLib' 'https://github.com/samtools/htslib.git'
	echo "Installing HTSLib..."
	rm -Rf htslib
	git clone https://github.com/samtools/htslib.git
	export HTSLIB_SRC=$PWD/htslib
	cd htslib
		make
	cd ..
}

install_libBigWig() {
	echo "Installing LibBigWig..."
	git clone https://github.com/dpryan79/libBigWig.git
	cd libBigWig
		make
	cd ..
}

#######################
## Brew
#######################
(uname -a | grep Darwin &> /dev/null) && [ ! -x "$(which brew)" ] && install_brew

#######################
## Zlib
#######################

echo "Checking for Zlib..."
echo '#include <zlib.h>' > tmp.c
echo 'int main(int c, char **v) {&deflate;}' >> tmp.c
(gcc tmp.c -o tmp.a &> tmp.o) || install_zlib

#######################
## Curl
#######################

echo "Checking for Curl..."
echo '#include <curl/multi.h>' > tmp.c
echo 'int main(int c, char **v) {&curl_multi_setopt;}' >> tmp.c
(gcc tmp.c -o tmp.a &> tmp.o) || install_zlib

#######################
## GSL
#######################

echo "Checking for GSL..."
echo "#include <gsl/gsl_sf_bessel.h>" > tmp.c
echo 'int main(int c, char **v) {&gsl_sf_bessel_J0;}' >> tmp.c
(gcc tmp.c -o tmp.a &> tmp.o) || install_gsl

#######################
## Git
#######################
echo "Checking for Git..."
if [ ! -x "$(which git)" ]; then
	install_git
fi

#######################
## HTSlib
#######################
echo "Checking for htslib"
if [ -z "$HTS_SRC" ];
then
	install_htslib
fi

#######################
## LibBigWig 
#######################
echo "Checking for LibBigWig src..."
if [ -z "$LIBBIGWIG_SRC" ];
then
	install_libBigWig
fi

#######################
## WiggleTools 
#######################
echo "Compiling WiggleTools"
make

rm tmp*
