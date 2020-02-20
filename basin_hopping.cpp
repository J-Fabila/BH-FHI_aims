/// De ser necesario agregar una celda como en VASP, que se ponga en geometry.in. En el 'else' de crystal==0 sustituye el '>' por '>>'
#include"atomicpp.h"
string Simbolo_1, Simbolo_2, file_name, command, aux,geometry_file, initialization_file, outputfile, i_str, E_str, tag;
int continue_alg,  Ncore, randomness, kick, iteraciones,swap_step, contenido, m, N_Simbolo_1, N_Simbolo_2, count, resto, crystal;
float step_width, Temperature, Energy, Energia, EnergiaAnterior, k_BT ;
Cluster clus_1, clus_2, clus;
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
cout<<"Continue = "<<continue_alg<<endl;
int i = 1;

if(continue_alg==1){
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //                                      RESTART ALGORITHM                                         //
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
      string iteration_counter_i ="cd ";
             iteration_counter_i+=file_name;
             iteration_counter_i+=" ; for i in $(ls coord*xyz ); do head -2 $i | tail -1 | awk '{ print $2 }' ; done | sort -n  | tail -1";
  i=int_pipe(iteration_counter_i,1);

      string iteration_counter_m ="cd ";
             iteration_counter_m+=file_name;
             iteration_counter_m+="/rejected ; ls coord*xyz 2> /dev/null | wc -l ";
  m=int_pipe(iteration_counter_m,0);
}
else{
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //                                        BEGIN ALGORITHM                                         //
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
   // Creates work directory
   command ="if [ -d "+file_name+" ] ; then mv "+file_name+" other_"+file_name;
   command+=" ; fi ; mkdir "+file_name+" ; cd "+file_name+"  ; mkdir rejected ;";
   //system(command.c_str()); command.clear(); command="cd input ; echo 'running command' >> run.sh "; system(command.c_str()); command.clear();
   command+=" cp ../input/* .";
   system(command.c_str());
   i=1; m=0;
   contenido=0;
   count=1;
   while(contenido!=1)
   {
      if(initialization_file.length() > 5 && count==1)
      {
      cout<<"Initialization file = "<<initialization_file<<endl;
      //Inicializa archivo geometry.in
        command ="cp ";
        command+=initialization_file;
        command+=" "+file_name+"/geometry.in";
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
            geometry_file=file_name+"/geometry.in";
            clus.print_fhi(geometry_file);
         }
         else{
            geometry_file.clear();
            geometry_file=file_name+"/geometry.tmp";
            clus.print_fhi(geometry_file);
            command.clear(); command="cd "+file_name+" ; cat crystal.in > geometry.in ; cat geometry.tmp >> geometry.in ; rm geometry.tmp ";
            system(command.c_str());
         }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      }
      command.clear();
      command="cd "+file_name+" ; ./run.sh";
      system(command.c_str());
      command.clear();
      command="cd "+file_name+" ; grep 'Have a nice day' output.out | wc -l";
      contenido=int_pipe(command.c_str());
   }
   command.clear();
   command="cd "+file_name+" ; grep \" | Total energy of the DFT \" output.out | cut -d \":\" -f 2 | cut -d \"e\" -f 1 ";
   Energy=float_pipe(command.c_str());
   i_str=to_string(i);
   E_str=to_string(Energy);
   command.clear();
   command=file_name+"/geometry.in.next_step";
   clus.read_fhi(command); command.clear(); command=file_name+"/coordinates1.xyz";
   tag.clear();
   tag=" Iteration "+i_str+" -----> Energy = "+E_str+" eV ";
   clus.print_xyz(command,tag);
   command.clear(); command="cd "+file_name+" ; mv output.out output1.out ; ";
   command+=" mv geometry.in geometry1.in ; echo Step Energy [eV] >> energies.txt ; echo 1  "+E_str+" >> energies.txt";
   system(command.c_str());

   cout<<" --> Relaxation of initial configuration: DONE! "<<endl<<endl;
   cout<<"=================================================================="<<endl;
   cout<<"BH-DFT routine starts here! "<<endl;
   cout<<"Note: "<<endl;
   cout<<"For monometallic clusters: only random xyz moves will be applied "<<endl;
   cout<<"For bimetallic clusters  : 1 atomic swap will be performed after "<<swap_step<<" moves "<<endl;
   cout<<"==================================================================="<<endl<<endl;
   i=2;
}
while(i+m <= iteraciones)
{
  resto=i%swap_step;
  geometry_file.clear();
  geometry_file=file_name+"/geometry.in.next_step";
  i_str.clear(); E_str.clear();
  i_str=to_string(i);
  // Get cluster coordinates from output file
  if(N_Simbolo_2>0) //es bimetalico
  {
    clus_1=extract(geometry_file,Simbolo_1);
    clus_2=extract(geometry_file,Simbolo_2);
    clus  =clus_1+clus_2;
  }
  else //es monometalico
  {
    clus=extract(geometry_file,Simbolo_1);
  }
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
    //aplica move;
    if(kick==0)
    {
      cout<<"Kicking without potential"<<endl;
      clus.kick(step_width);
    }
    else
    {
      cout<<"Kicking with lennard potential"<<endl;
      clus.kick_lennard(step_width);
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(crystal==0){
      geometry_file.clear();
      geometry_file=file_name+"/geometry.in";
      clus.print_fhi(geometry_file);
  }
  else{
      geometry_file.clear();
      geometry_file=file_name+"/geometry.tmp";
      clus.print_fhi(geometry_file);
      command.clear(); command="cd "+file_name+" ; cat crystal.in > geometry.in ; cat geometry.tmp >> geometry.in ; rm geometry.tmp ";
      system(command.c_str());
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  command.clear();
  command="cd "+file_name+" ; ./run.sh";
  system(command.c_str());
  command.clear();
  command="cd "+file_name+" ; grep 'Have a nice day' output.out | wc -l";
  contenido=int_pipe(command);
// Starting randomly if last configuration fails
  while (contenido!=1)
  {
  cout<<" --> SCF failed. Starting again from randomly generated structure! "<<endl;
// sera necesario borrar geometry.in???
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
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(crystal==0){
      geometry_file.clear();
      geometry_file=file_name+"/geometry.in";
      clus.print_fhi(geometry_file);
  }
  else{
      geometry_file.clear();
      geometry_file=file_name+"/geometry.tmp";
      clus.print_fhi(geometry_file);
      command.clear(); command="cd "+file_name+" ; cat crystal.in > geometry.in ; cat geometry.tmp >> geometry.in ; rm geometry.tmp ";
      system(command.c_str());
  }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  command.clear();
  command="cd "+file_name+" ;  ./run.sh";
  system(command.c_str());
  command.clear();
  command="cd "+file_name+" ; grep 'Have a nice day' output.out | wc -l";
  contenido=int_pipe(command.c_str());
  }
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//                                         SAVE ENERGIES                                          //
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
EnergiaAnterior=Energy;
command.clear();
command="cd "+file_name+" ; grep \" | Total energy of the DFT \" output.out | cut -d \":\" -f 2 | cut -d \"e\" -f 1 ";
Energy=float_pipe(command.c_str());
E_str=to_string(Energy);
command.clear();
command=file_name+"/geometry.in.next_step";
clus.read_fhi(command); command.clear();
command=file_name+"/coordinates"+i_str+".xyz";
tag.clear(); tag=" Iteration "+i_str+" -----> Energy = "+E_str+" eV ";
clus.print_xyz(command,tag);
command.clear(); command="mv "+file_name+"/output.out "+file_name+"/output"+i_str+".out";
system(command.c_str());
command.clear(); command="mv "+file_name+"/geometry.in "+file_name+"/geometry"+i_str+".in";
system(command.c_str());
command.clear(); command="echo "+i_str+"  "+E_str+" >> "+file_name+"/energies.txt";
system(command.c_str() );
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//                                     Metropolis Monte-Carlo                                     //
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
k_BT = 0.00008617 * Temperature;
if (pow(2.71,(EnergiaAnterior-Energy)/k_BT) > random_number(0,1))
{
  cout<<"--> Basin Hopping MC criteria: Energy accepted! "<<endl;
  cout<<"--> Finished iteration "<<i<<endl;
  command.clear(); command="tail -"+i_str+" "+file_name+"/energies.txt |  sort -nk2 > "+file_name+"/sorted.txt";
  system(command.c_str());
  i++;
}
else
{
  cout<<"--> Basin Hopping MC criteria: Energy rejected!"<<endl;
  command.clear(); command="mv "+file_name+"/output"+i_str+".out "+file_name+"/rejected/outputrejected"+i_str+".out";
  system(command.c_str());
  command.clear(); command="mv "+file_name+"/geometry"+i_str+".in "+file_name+"/rejected/geometryrejected"+i_str+".in";
  system(command.c_str());
  m++;
}
}
return 0;
}
