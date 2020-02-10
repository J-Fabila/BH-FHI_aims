#Download files
echo "
"
echo -e "\e[34m Downloading files ... \e[0m"
echo "
"
sleep 1.5
wget https://raw.githubusercontent.com/J-Fabila/Atomic.h/master/atomicpp.h
wget https://raw.githubusercontent.com/J-Fabila/BH-FHI_aims/master/basin_hopping.cpp
wget https://raw.githubusercontent.com/J-Fabila/BH-FHI_aims/master/input.bh
wget https://raw.githubusercontent.com/J-Fabila/BH-VASP/master/Instruction_Manual.pdf
g++ -o basin_hopping basin_hopping.cpp
echo -e "\e[34m Compiled ... \e[0m"
echo "
"
sleep 1.7

echo "
"
echo -e "\e[34m Creating directories ... \e[0m"
echo "
"
sleep 1.5
mkdir input
cd input
wget https://raw.githubusercontent.com/J-Fabila/Atomic.h/master/input/control.in
wget https://raw.githubusercontent.com/J-Fabila/Atomic.h/master/input/geometry.in
cd ..
rm $0
echo -e "
\e[34m                     Done!
\e[0m"
sleep 1.8
clear
