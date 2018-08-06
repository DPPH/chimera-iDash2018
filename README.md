# chimera-iDash2018

installation prerequisites (ubuntu 18.04)
```sh
sudo apt-get install libntl-dev libgmp-dev cmake
```


compile instructions
```sh
git submodule init
git submodule update

#then use Clion editor, or compile manually as below
mkdir build
cd build; 
cmake ../fhe
make
```
