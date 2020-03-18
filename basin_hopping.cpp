#include<stdio.h>
#include<stdio.h>
#include<unistd.h>
#include"atomicpp.h"
string Simbolo_1, Simbolo_2, file_name, command, aux,geometry_file, initialization_file, outputfile, i_str, E_str, tag;
int continue_alg,  Ncore, randomness, kick, iteraciones,swap_step, contenido, m, N_Simbolo_1, N_Simbolo_2, count, fail_counter=0, resto, failed_max,crystal;
float step_width, Temperature, Energy, Energia, EnergiaAnterior, k_BT, damp ;
float x_min,y_min,z_min,x_max,y_max,z_max;
Cluster clus_1, clus_2, clus, c_aux;
Crystal cristal;
float dist=0.0;
int main(int argc, char *argv[]){
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//                                    Gets data from input.bh                                     //
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
Simbolo_1=string_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 2 | cut -d ':' -f 1 ");
Simbolo_2=string_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 3 | cut -d ':' -f 1 ");
N_Simbolo_1=int_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 2 | cut -d ':' -f 2 | cut -d ']' -f 1 ");
N_Simbolo_2=int_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 3 | cut -d ':' -f 2 | cut -d ']' -f 1 ");
continue_alg=int_pipe("grep 'continue' input.bh | awk '{print $3}' ");
initialization_file=string_pipe("grep 'initialization_file' input.bh | awk '{print $3}' ");
randomness=int_pipe("grep 'randomness' input.bh | awk '{print $3}' ");
kick=int_pipe("grep 'kick_type' input.bh | awk '{print $3}' ");
file_name=string_pipe("grep 'file_name' input.bh | awk '{print $3}' ");
step_width=float_pipe("grep 'step_width' input.bh | awk '{print $3}' "); //o double?
Temperature=float_pipe("grep 'temperature_K' input.bh | awk '{ print $3 }' "); //o double?
Ncore=int_pipe("grep 'Ncore' input.bh | head -1 | awk '{print $3}' ");
iteraciones=int_pipe("grep 'iterations' input.bh | awk '{ print $3 }' ");
swap_step=int_pipe("grep 'swap_step' input.bh | awk '{ print $3 }' ");
crystal=int_pipe("cd input ; if [ -f crystal.in ]  ; then echo 1  ;  fi ");
failed_max=3;
damp=0.0;
cout<<Simbolo_1<<"  "<<Simbolo_2<<"  "<<kick<<" "<<randomness<<endl;
//Put if crystal==1 and a grep in input.bh were is indicated supported or not??
if(crystal==1)  //Esto sustituye tener que poner [x_min,x_max]; [y_min,y_max]... en el input
{
  cristal.read_fhi("input/crystal.in");
  x_min=cristal.x_min();
  x_max=cristal.x_max();
  y_min=cristal.y_min();
  y_max=cristal.y_max();
  z_min=cristal.z_min();
  z_max=cristal.z_max();
}
int i = 1;

if(continue_alg==1)
{
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //                                      RESTART ALGORITHM                                         //
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
      chdir(file_name.c_str());
      string iteration_counter_i ="ls coord*xyz | wc -l";
      i=int_pipe(iteration_counter_i);

      string iteration_counter_m ="ls ./rejected/coord*xyz | wc -l ";
      m=int_pipe(iteration_counter_m);
}


else
{
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //                                        BEGIN ALGORITHM                                         //
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
   // Creates work directory
   command ="if [ -d "+file_name+" ] ; then mv "+file_name+" other_"+file_name; ///// Modify to be in c++ code??
   command+=" ; fi ; mkdir "+file_name+" ; cd "+file_name+"  ; mkdir rejected ;";
   command+=" cp ../input/* .";
   system(command.c_str());
   command.clear();
   chdir(file_name.c_str());
   i=1; m=0;
   contenido=0;
   count=1;
   while(contenido!=1)
   {
      if(initialization_file.length() > 5 && count==1)
      {
      //Inicializa archivo geometry.in
        command ="cp ";
        command+=initialization_file;
        command+="geometry.in";
        system(command.c_str());
        command.clear();
        count++;
      }
      else
      {
      //Empieza desde 0
         if(N_Simbolo_2>0)  //Esto quiere decir que es bimetalico
         {
            if(randomness==1)  // fully random
            {
               clus.srand_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
            }
            else //pseudoaleatorio
            {
               clus.rand_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
            }
         }
         else //es monometalico
         {
            if(randomness==1)  // fully random
            {
               clus.srand_generator(Simbolo_1,N_Simbolo_1);
            }
            else //pseudoaleatorio
            {
               clus.rand_generator(Simbolo_1,N_Simbolo_1);
            }
         }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         if(crystal==0){
            geometry_file.clear();
            geometry_file="geometry.in";
            clus.print_fhi(geometry_file);
         }
         else{
            clus.centroid();
            clus.move((x_max-x_min)/2.0+random_number(-dist,dist),(y_max-y_min)/2.0+random_number(-dist,dist),z_max-clus.z_min());
            geometry_file.clear();
            geometry_file="geometry.tmp";
            clus.print_fhi(geometry_file);
            command.clear();
            command="cat crystal.in > geometry.in ; cat geometry.tmp >> geometry.in ; rm geometry.tmp ";
            system(command.c_str());
            command.clear();
         }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      }

      // RUN FIRST CONFIGURATION:
      command="./run.sh";
      system(command.c_str());
      command.clear();
      command="grep 'Have a nice day' output.out | wc -l";
      contenido=int_pipe(command.c_str());
      command.clear();
   }


   /// Store energy and optimized geometry:
   command="grep \" | Total energy of the DFT \" output.out | awk '{print $12}' ";
   Energy=double_pipe(command.c_str());
   i_str=to_string(i);
   E_str=string_pipe(command); //Better for Energies with all the value
   command.clear();

   // get initial, relaxed geometry and store it as xyz (first xyz!)
   command="geometry.in.next_step";
   clus.read_fhi(command);
   command.clear();
   command="coordinates1.xyz";
   tag.clear();
   tag=" Iteration "+i_str+" -----> Energy = "+E_str+" eV ";
   clus.print_xyz(command,tag);
   command.clear();


   // store first output and geometry files:
   command="mv output.out output1.out";
   system(command.c_str());
   command.clear();

   command="mv geometry.in geometry1.in";
   system(command.c_str());
   command.clear();

   command="echo 'Step ----> Energy[eV]' >> energies.txt ; echo '1 ---->' "+E_str+" >> energies.txt";
   system(command.c_str());
   command.clear();



   // NOW THE BASIN HOPPING WILL BEGIN:

   cout<<" --> Relaxation of initial configuration: DONE! "<<endl;
   cout<<"================================================================================================"<<endl;
   cout<<"BH-DFT routine starts here! "<<endl;
   cout<<"Note: "<<endl;
   cout<<"For monometallic clusters: only random xyz moves will be applied "<<endl;
   cout<<"For bimetallic clusters  : 1 atomic swap will be performed after "<<swap_step<<" moves "<<endl;
   cout<<"================================================================================================"<<endl;


   // set up BH-DFT iteration counter to 2:
   i=2;


// finished else for initial, randomly generated structure!
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
// BASIN HOPPING LOOP:
// get coordinates and apply kick
////////////////////////////////////////////////////////////////////////////////////////////////////////



while(i+m <= iteraciones)
{
  // when to do atomic swap:
  resto=i%swap_step;

  // cleaning memory
  i_str.clear();
  E_str.clear();
  i_str=to_string(i);


  // For bimetallics:
  if(N_Simbolo_2>0)
  {
   // get relaxed coordinates from initial configuration (i=1):
   geometry_file.clear();
   geometry_file="geometry.in.next_step";
   command.clear();
   command="cat "+geometry_file;
   system(command.c_str());
   clus_1=extract(geometry_file,Simbolo_1);
   clus_2=extract(geometry_file,Simbolo_2);

   clus  =clus_1+clus_2;
  }


  // for monometallics:
  else
  {
   // get relaxed coordinates from initial configuration (i=1):
   geometry_file.clear();
   geometry_file="geometry.in.next_step";
   command.clear();
   command="cat "+geometry_file;
   system(command.c_str());
   clus=extract(geometry_file,Simbolo_1);
  }




  /////////////////////////////////////
  // MOVE CLASSES: KICK VS ATOMIC SWAP
  /////////////////////////////////////


  // Applies swap for bimetallic case by calling "atomic.cpp"
  if(resto==0 && N_Simbolo_2 > 0)
  {
   cout<<"Applying atomic swap"<<endl;
   clus.type = "bimetallic";

    if(N_Simbolo_1>=N_Simbolo_2)
    {
       clus.swap(N_Simbolo_2);
    }
    else
    {
       clus.swap(N_Simbolo_1);
    }
  }


  // Applies kick for monometallic cases and non-swapping ones
  // by calling "atomic.cpp"
  else
  {
    if(kick==0)
    {
      // random kick:
      cout<<"Kicking without potential"<<endl;
      clus.kick(step_width);
    }


    else
    {
      // damped kick (jorge):
      cout<<"Kicking with lennard potential"<<endl;
      clus.kick_lennard(step_width); //Need to put range of -1 to 1 of the random value * step width
    }
  }



  // Update geometries if in the gas-phase:
  if(crystal==0)
   {
      geometry_file.clear();
      geometry_file="geometry.in";
      clus.print_fhi(geometry_file);
   }

  //.. or periodic case:
  else
   {
      clus.centroid();
      clus.move((x_max-x_min)/2.0+random_number(-dist,dist),(y_max-y_min)/2.0+random_number(-dist,dist),z_max-clus.z_min() );
      geometry_file.clear();
      geometry_file="geometry.tmp";
      clus.print_fhi(geometry_file);
      command="cat crystal.in > geometry.in ; cat geometry.tmp >> geometry.in ; rm geometry.tmp ";
      system(command.c_str());
      command.clear();
   }


  ////////////////////////////////////////////
  // RUNNING FHI-AIMS SECOND CONFIGURATION....:
  command="./run.sh";
  system(command.c_str());
  command.clear();

  // Checking that geometry has converged:
  command="grep 'Have a nice day' output.out | wc -l";
  contenido=int_pipe(command);
  command.clear();

  //***********************************************************************************
  // If structure DID NOT converged, then:
  while (contenido!=1)
  {
  cout<<" --> SCF failed.  "<<endl;


  if(fail_counter<failed_max)
  {
     cout<<" Kicking ... again"<<endl;
     fail_counter++;
     geometry_file.clear();
     geometry_file="coordinates"+i_str+".xyz";
     c_aux.read_xyz(geometry_file);
     geometry_file.clear();
     geometry_file="aux.fhi";
     c_aux.print_fhi(geometry_file);

     // for bimetallics:
     if(N_Simbolo_2>0)
     {
       clus_1=extract(geometry_file,Simbolo_1);
       clus_2=extract(geometry_file,Simbolo_2);
       clus  =clus_1+clus_2;
     }

     // for monometallics:
     else
     {
       clus=extract(geometry_file,Simbolo_1);
     }


     system("rm aux.fhi");


     // Applies swap or kick
     if(resto==0 && N_Simbolo_2 > 0)
     {
      cout<<"Applying swap step"<<endl;
      clus.type = "bimetallic";
       if(N_Simbolo_1>=N_Simbolo_2)
       {
          clus.swap(N_Simbolo_1);
       }
       else
       {
          clus.swap(N_Simbolo_2);
       }
     }


     else
     {
       if(kick==0)
       { // increase kick with damp to explore PES
         cout<<"Kicking without potential"<<endl;
         clus.kick(step_width+(fail_counter*damp));
       }


       else
       { // increase kick with damp to explore PES
         cout<<"Kicking with lennard potential"<<endl;
         clus.kick_lennard(step_width+(fail_counter*damp));
       }
     }
  }



  else
  {
     cout<<" Starting again from randomly generated structure!"<<endl;
     if(N_Simbolo_2>0)
     {
         if(randomness==1)  // fully random
         {
            clus.srand_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
         }
         else //pseudoaleatorio
         {
            clus.rand_generator(Simbolo_1,N_Simbolo_1,Simbolo_2,N_Simbolo_2);
         }
     }
     else //es monometalico
     {
         if(randomness==1)  // fully random
         {
            clus.srand_generator(Simbolo_1,N_Simbolo_1);
         }
         else //pseudoaleatorio
         {
            clus.rand_generator(Simbolo_1,N_Simbolo_1);
         }
     }
  }


  if(crystal==0){
      geometry_file.clear();
      geometry_file="geometry.in";
      clus.print_fhi(geometry_file);
  }

  else{
      clus.centroid();
      clus.move((x_max-x_min)/2.0+random_number(-dist,dist),(y_max-y_min)/2.0+random_number(-dist,dist),z_max-clus.z_min() );
      geometry_file.clear();
      geometry_file="geometry.tmp";
      clus.print_fhi(geometry_file);
      command="cat crystal.in > geometry.in ; cat geometry.tmp >> geometry.in ; rm geometry.tmp ";
      system(command.c_str());
      command.clear();
  }


  // THIRD RUN FOR NON-CONVERGED STRUCTURES

  command="./run.sh";
  system(command.c_str());
  command.clear();

  // IF CONVERGED....
  command="grep 'Have a nice day' output.out | wc -l";
  contenido=int_pipe(command.c_str());
  command.clear();
  }

//***********************************************************************************************************************************




//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//                                         SAVE ENERGIES                                          //
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

EnergiaAnterior=Energy;

command="grep \" | Total energy of the DFT \" output.out | awk '{print $12}' ";
Energy=float_pipe(command.c_str());
E_str=string_pipe(command);
command.clear();

/////Perhaps this part can be inside Metropolis Monte-Carlo (in step of mv outi mv out)
command="geometry.in.next_step";
clus.read_fhi(command);
command.clear();


command="coordinates"+i_str+".xyz";
tag.clear();
tag=" Iteration "+i_str+" -----> Energy = "+E_str+" eV ";
clus.print_xyz(command,tag);
command.clear();


command="mv output.out output"+i_str+".out";
system(command.c_str());
command.clear();


command="mv geometry.in geometry"+i_str+".in";
system(command.c_str());
command.clear();



//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//                                     Metropolis Monte-Carlo                                     //
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
k_BT = 0.00008617 * Temperature;
if (pow(2.71,(EnergiaAnterior-Energy)/k_BT) > random_number(0,1))
{
  cout<<"--> Basin Hopping MC criteria: Energy accepted! "<<endl;
  cout<<"--> Finished iteration "<<i<<endl;
  command="echo "+i_str+" '---->'  "+E_str+" >> energies.txt";
  system(command.c_str());
  command.clear();
  command="tail -"+i_str+" energies.txt |  sort -nk3 > sorted.txt";
  system(command.c_str());
  command.clear();
  i++;
  fail_counter=0;
}


else
{
  cout<<"--> Basin Hopping MC criteria: Energy rejected!"<<endl;
  m++;
  string m_str = to_string(m);
  command="mv output"+i_str+".out ./rejected/outputrejected"+m_str+".out";
  system(command.c_str());
  command.clear();
  command="mv geometry"+i_str+".in ./rejected/geometryrejected"+m_str+".in";
  system(command.c_str());
  command.clear();
  m_str.clear();
  fail_counter=0;
}


} // END OF BH-LOOP
return 0;
}
