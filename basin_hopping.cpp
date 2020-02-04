#include<atomicpp.h>

int continue_alg, randomness, kick, iteraciones,swap_step, contenido;
float step_width, temperature ;
string file_name, command;
int main(int argc, char *argv[]){

//################################################################################################
//#                                    Gets data from input.bh                                   #
//################################################################################################

continue_alg = stoi(argv[1]);
randomness   = stoi(argv[2]);
kick         = stoi(argv[3]);
iteraciones  = stoi(argv[4]);
step_width   = stof(argv[5]);
temperature  = stof(argv[6]);
file_name    = argv[7]);
swap_step    = stoi(argv[8]);
sel          = stoi(argv[9]);
Nat          = stoi(argv[10]);
n            = stoi(argv[11]);
Ncore        = stoi(argv[12]);
Simbolo_1    = argv[12];
N_Simbolo_1  = stoi(argv[13]);
path         = argv[14];

if(argc == 17 )
{
   Simbolo_2=argv[14];
   N_Simbolo_2=stoi(argv[15]);
}

//################################# CREATES Work Directory #####################################
command="if [ -d $file_name ] ; then mv $file_name other_$file_name ; fi ; cp -r input $file_name";


//##############################################################################################
//#                                         BEGIN ALGORITHM                                    #
//##############################################################################################


system("mkdir rejected");
contenido=0;
count=1;
while(contenido!=1)
{

  if(initialization_file==1 && count==1)
  {
    //Inicializa archivo geometry.in
    count++;
  }
  else
  {
    //Empieza desde 0
    if(n>3)  //Esto quiere decir que es bimetalico
    {
      if(randomness==1)  // fully random
      {
        // clus.SRandomGenerator(...)
      }
      else //pseudoaleatorio
      {
        // clus.RandomGenerator(...)
      }
    }
    else //es monometalico
    {
        if(randomness==1)  // fully random
        {
          // clus.SRandomGenerator(...)
        }
        else //pseudoaleatorio
        {
          // clus.RandomGenerator(...)
        }
    }
  }
  if(sel==1)
  {
    //algoritmo que haga lo propio cuando hay selective dynamics
  }
  system("./run.sh");
  system("grep 'reached required ' OUTCAR | wc -l > contenido.txt");
  ifstream cont...
  //read contenido.txt
}

/////////  MODIFICAR ESTO PARA FHI AIMS (AHORITA ESTA PARA vasp)
Energia=$(tail -1 OSZICAR | awk '{print $5 }')                 #Extrae  la energia del OSZICAR
echo "1     $Energia " >> CONTCAR1      #Escribe la energia y num de iteracion en CONTCAR1
N=$(wc -l CONTCAR | awk '{ print $1 }' )         #Cuenta el numero de lineas que tiene CONTCAR
tail -$(($N-1)) CONTCAR >> CONTCAR1       #Mueve la informacion a CONTCAR1 con el nuevo titulo
mv POSCAR POSCAR1                                                      #Mueve POSCAR a POSCAR1
rm CHG CHGCAR DOSCAR EIGENVAL XDATCAR IBZKPT OSZICAR PCDAT REPORT WAVECAR *.xml CONTCAR
echo "Step  Energy" >> energies.txt;
echo "1     $Energia" >> energies.txt

cout<<" --> Relaxation of initial configuration: DONE! "<<endl<<endl;
cout<<"=================================================================="<<endl;
cout<<"BH-DFT routine starts here! "<<endl;
cout<<"Note: "<<endl;
cout<<"For monometallic clusters: only random xyz moves will be applied "<<endl;
cout<<"For bimetallic clusters  : 1 atomic swap will be performed after 10 moves "<<endl;
cout<<"==================================================================="<<endl<<endl;

i=2

while(i+m < iteraciones)
{
  resto=i%swap_step;
  //extract cluster
  if(n>3) //es bimetalico
  {
    cluster=extract(Simbolo_1,outputfile);
    cluster=extract(Simbolo_2,outputfile);
    cluster=cluster+cluster;
  }
  else //es monometalico
  {
    cluster=extract(Simbolo_1,outputfile);
  }
  if(resto==0)
  {
    //aplica swap;
  }
  else
  {
    //aplica move;
    if(kick==0)
    {
      //aplica tipo depatada =0
    }
    else
    {
      //aplica la otra patada
    }
  }
  system("./run.sh");
  system("grep 'reached required ' OUTCAR | wc -l > contenido.txt");
  ifstream cont...
  //read contenidoi=contenido.txt

  while (contenido!=1)
  {
  cout<<" --> SCF failed. Starting again from randomly generated structure! "<<endl;

  if(n>3)  //Esto quiere decir que es bimetalico
  {
    if(randomness==1)  // fully random
    {
      // clus.SRandomGenerator(...)
    }
    else //pseudoaleatorio
    {
      // clus.RandomGenerator(...)
    }
  }
  else //es monometalico
  {
      if(randomness==1)  // fully random
      {
        // clus.SRandomGenerator(...)
      }
      else //pseudoaleatorio
      {
        // clus.RandomGenerator(...)
      }
  }

  system("./run.sh");
  system("grep 'reached required ' OUTCAR | wc -l > contenido.txt");
  ifstream cont...
  //read contenidoi=contenido.txt
  }

  ################################################################################################
  #                                     Save  configuration                                      #
  ################################################################################################

     EnergiaAnterior=$(echo $Energia)
     Energia=$(tail -1 OSZICAR | awk '{print $5 }')    #Extrae energia del OSZICAR

     echo "$i     $Energia " >> CONTCAR$i
     N=$(wc -l CONTCAR | awk '{ print $1 }' )
     tail -$(($N-1)) CONTCAR >> CONTCAR$i  #Renombra los archivos y les agrega la energia
     mv POSCAR POSCAR$i
     rm CHG CHGCAR DOSCAR OSZICAR EIGENVAL XDATCAR IBZKPT  PCDAT REPORT WAVECAR *.xml CONTCAR

  ################################################################################################
  #                                  Metropolis Monte-Carlo                                      #
  ################################################################################################

     accepted=$(python ../programs/metropolis.py $EnergiaAnterior $Energia $Temperature)

        if [ $accepted = true ]
        then   # La energía ha sido aceptada

          echo "--> Basin Hopping MC criteria: Energy accepted! "
          echo "--> Finished iteration $i"
          head -1 CONTCAR$i >> energies.txt
          tail -$i energies.txt |  sort -nk2 > sorted.txt
          i=$(($i+1)) #Convergencia lograda, aumenta en 1 el contador

        else

          echo "--> Basin Hopping MC criteria: Energy rejected!"
          mv CONTCAR$i rejected/CONTCARrejected$i
          mv POSCAR$i rejected/POSCARrejected$i
          m=$(($m+1))

        fi
}
return 0;
}
