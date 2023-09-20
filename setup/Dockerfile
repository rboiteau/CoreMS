FROM debian:bullseye

RUN apt update && apt upgrade -y
RUN apt install clang build-essential pkg-config autoconf libtool python3-dev libhdf5-dev libnetcdf-dev -y
RUN apt install wget -y
RUN apt install zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libsqlite3-dev libreadline-dev libffi-dev curl libbz2-dev -y
RUN apt install gdb lldb vim -y
RUN apt install git -y


RUN apt install apt-transport-https dirmngr gnupg ca-certificates -y
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF

RUN echo "deb https://download.mono-project.com/repo/debian stable-buster main" | tee /etc/apt/sources.list.d/mono-official-stable.list

RUN apt update && apt install mono-complete nuget -y

RUN wget https://dot.net/v1/dotnet-install.sh -O dotnet-install.sh
RUN chmod +x dotnet-install.sh
RUN ./dotnet-install.sh --version latest
RUN apt update

WORKDIR /

RUN git clone https://github.com/deweycw/CoreMS.git
WORKDIR /CoreMS
RUN apt update
RUN apt install python3-pip -y

RUN git checkout stable-docker
RUN python3 -m pip install -U pip
RUN python3 -m pip install wheel
RUN python3 -m pip install numpy
RUN python3 -m pip install netCDF4==1.5.3
RUN python3 -m pip install -r requirements-1.txt
RUN python3 ./setup.py install


