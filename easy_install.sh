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

install_mysql () {
	echo "Installing mysql...."
	if [ -x "$(which apt-get)" ]; then
		apt-get install mysql-server mysql-client
	elif [ -x "$(which brew)" ]; then 
		brew install mysql
	else
		if [ ! -x "$(which cmake)" ];
		then
			echo "Installing cmake..."
			rm -Rf cmake*
			confirmed_download 'cmake' 'http://www.cmake.org/files/v3.0/cmake-3.0.0.tar.gz'
			cd cmake-3.0.0
				./configure
				make
				sudo make install
			cd ..
			export PATH=$PATH
		fi 
		rm -Rf mysql*
		confirmed_download 'mysql' 'http://sourceforge.net/projects/mysql.mirror/files/MySQL%205.6.19/mysql-5.6.19.tar.gz/download'
		cd mysql*
			cmake .
			make 
			sudo make install
		cd ..
	fi
}

install_ssl () {
	if [ -x "$(which apt-get)" ]; then
		apt-get install openssl
	elif [ -x "$(which brew)" ]; then
		brew install openssl
	else 
		rm -Rf openssl*
		confirmed_download 'OpenSSl' 'ftp://ftp.openssl.org/source/openssl-1.0.1h.tar.gz'
		cd openssl*
			./config
			make
			sudo make install
		cd ..
	fi
}

install_png () {
	if [ -x "$(which apt-get)" ]; then
		apt-get install libpng-dev
	elif [ -x "$(which brew)" ]; then
		brew install libpng
	else 
		rm -Rf libpng*
		confirmed_download 'libpng' 'http://prdownloads.sourceforge.net/libpng/libpng-1.6.12.tar.gz?download'
		cd libpng
			./configure
			make
			sudo make install
		cd ..
	fi
}

install_gsl() {
	if [ -x "$(which apt-get)" ]; then
		apt-get install libgsl0 libgsl0-dev
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

install_tabix() {
	prompt_download 'Tabix' 'https://github.com/samtools/tabix.git'
	echo "Installing Tabix..."
	rm -Rf tabix
	git clone https://github.com/samtools/tabix.git
	export TABIX_SRC=$PWD/tabix
	cd tabix
		make
	cd ..
}

install_kent() {
	prompt_download 'Kent userApps' 'genome-source.cse.ucsc.edu/kent.git'
	echo "Installing Kent source..."
	rm -Rf userApps
	git archive --format=zip -9 --remote=git://genome-source.cse.ucsc.edu/kent.git beta src/userApps > userApps.zip
	unzip -d userApps -j userApps.zip
	rm userApps.zip
	cd userApps
		make fetchSource
		make libs
		export KENT_SRC=$PWD/kent/src
		mkdir bin
		cd kent/src/utils/bigWigCat
			make
		cd ../../../..
	cd ..
}

#######################
## Brew
#######################
(uname -a | grep Darwin &> /dev/null) && [ ! -x "$(which brew)" ] && install_brew

#######################
## png
#######################

echo "Checking for OpenSSL..."
echo '#include <openssl/ssl.h>' > tmp.c
echo 'int main(int c, char **v) {&SSL_CIPHER_description;}' >> tmp.c
(gcc tmp.c -o tmp.a &> tmp.o) || install_ssl

#######################
## png
#######################

echo "Checking for libpng..."
echo '#include <png.h>' > tmp.c
echo 'int main(int c, char **v) {&png_sig_cmp;}' >> tmp.c
(gcc tmp.c -o tmp.a &> tmp.o) || install_png

#######################
## MySQL
#######################

echo "Checking for mysql..."
echo "#include <mysql.h>" > tmp.c
echo 'int main(int c, char **v) {&mysql_connect;}' >> tmp.c
(gcc tmp.c -o tmp.a &> tmp.o) || install_mysql

#######################
## GSL
#######################

echo "Checking for GSL..."
echo "#include <gsl/gsl_sf_bessel.h>" > tmp
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
## Tabix
#######################
echo "Checking for tabix..."
if [ -z "$TABIX_SRC" ];
then
	install_tabix
fi

#######################
## Kent
#######################
echo "Checking for Kent src..."
if [ -z "$KENT_SRC" ];
then
	install_kent
fi

#######################
## WiggleTools 
#######################
echo "Compiling WiggleTools"
make

rm tmp*
