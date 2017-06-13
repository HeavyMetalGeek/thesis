# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
OF_USER_DIR=/user/f325j839/OpenFOAM/f325j839-3.0.0
FOAM_RUN=$OF_USER_DIR/run
FOAM_APPS=$OF_USER_DIR/applications
FOAM_BIN="/user/f325j839/OpenFOAM/OpenFOAM-3.0.0/bin"
FOAM_USER_BIN=$OF_USER_DIR/bin

SRC_OF2(){
    echo "source /apps/OpenFOAM/OpenFOAM-2.3.x/etc/bashrc"
    source '/apps/OpenFOAM/OpenFOAM-2.3.x/etc/bashrc'
}

SRC_OF3(){
    echo "source /apps/OpenFOAM/OpenFOAM-3.0.0/etc/bashrc"
    source '/apps/OpenFOAM/OpenFOAM-3.0.0/etc/bashrc'
}

PATH=$PATH:/apps/ParaView/ParaView-4.3.1-Linux-64bit/bin
PATH=$PATH:/user/f325j839/OpenFOAM/OpenFOAM-3.0.0/bin
