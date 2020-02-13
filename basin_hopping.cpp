#include"atomicpp.h"
string Simbolo_1, Simbolo_2;
int N_Simbolo_1, N_Simbolo_2, count, resto;
string initialization_file, outputfile;
int continue_alg,  Ncore, randomness, kick, iteraciones,swap_step, contenido, m;
float step_width, Temperature, Energy, Energia, EnergiaAnterior, k_BT ;
string file_name, command, aux;
Cluster clus_1, clus_2, clus;
int main(int argc, char *argv[])
{
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//                                    Gets data from input.bh                                     //
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
Simbolo_1=string_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 2 | cut -d ':' -f 1 ");
Simbolo_2=string_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 3 | cut -d ':' -f 1 ");
N_Simbolo_1=int_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 2 | cut -d ':' -f 2 | cut -d ']' -f 1 ");
N_Simbolo_2=int_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 3 | cut -d ':' -f 2 | cut -d ']' -f 1 ");
continue_alg=int_pipe("grep 'continue' input.bh | cut -d ' ' -f 3 ");
initialization_file=string_pipe("grep 'initialization_file' input.bh | awk '{print $3}' ");
randomness=int_pipe("grep 'randomness' input.bh | awk '{print $3}' ");
kick=int_pipe("grep 'kick_type' input.bh | awk '{print $3}' ");
file_name=string_pipe("grep 'file_name' input.bh | awk '{print $3}' ");
step_width=float_pipe("grep 'step_width' input.bh | awk '{print $3}' "); //o double?
Temperature=float_pipe("grep 'temperature_K' input.bh | awk '{ print $3 }' "); //o double?
Ncore=int_pipe("grep 'Ncore' input.bh | head -1 | awk '{print $3}' ");
iteraciones=system("grep 'iterations' input.bh | awk '{ print $3 }' ");
swap_step=system("grep 'swap_step' input.bh | awk '{ print $3 }' ");

if(continue_alg==1)//continue==1 significa que retoma un calculo anterior
{
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //                                      RESTART ALGORITHM                                         //
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

      string iteration_counter_i ="cd ";
             iteration_counter_i+=file_name;
             iteration_counter_i+=" ; for i in $(ls coord*xyz ); do head -1 $i | awk '{ print $1 }' ; done | sort -n  | tail -1";
  i=int_pipe(iteration_counter_i,1);

      string iteration_counter_m ="cd ";
             iteration_counter_m+=file_name;
             iteration_counter_m+="/rejected ; ls coord*xyz 2> /dev/null | wc -l ";
  m=int_pipe(iteration_counter_m,0);
}
else  //quiere decir que empieza desde scratch
{
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //                                        BEGIN ALGORITHM                                         //
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

   // Creates work directory
/*   command ="if [ -d ";
   command+=file_name;
   command+="] ; then mv ";
   command+=file_name;
   command+=" other_$";
   command+=file_name;
   command+=" ; fi ; cp -r input $";
   command+=file_name;
cout<<command<<endl;*/
command="mkdir "+file_name+" ; cd "+file_name+" ; mkdir rejected ";
   system(command.c_str());
   command.clear();
   command="cp input/* "+file_name;
   system(command.c_str());
   i=1; m=0;
   contenido=0;
   count=1;
   while(contenido!=1)
   {
      if(initialization_file>"alse" && count==1)
      {
      //Inicializa archivo geometry.in
        command ="cp ";
        command+=initialization_file;
        command+=" geometry.in";
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
         clus.print_fhi("geometry.tmp");
         command.clear();  command="cd "+file_name+"  cat ../geometry.tmp > geometry.in ; rm ../geometry.tmp";
system(command.c_str());
      }
      command.clear();
      command="cd "+file_name+"  ./run.sh";
      system(command.c_str());
      contenido=int_pipe("grep 'Have a nice day' output.out | wc -l");
   }
   Energy=double_pipe("grep \" | Total energy of the DFT / Hartree-Fock s.c.f. calculation \" output.out | cut -d \":\" -f 2 | cut -d \"e\" -f 1 ");
   system("mv geometry.in.next_step coordinates1.xyz");
   system("mv output.out output1.out");
   system("mv geometry.in geometry1.in");
   system("echo Step Energy [eV] >> energies.txt");
          command ="echo 1  ";
          command+=to_string(Energy);
          command+=" >> energies.txt";
   system(command.c_str());

   cout<<" --> Relaxation of initial configuration: DONE! "<<endl<<endl;
   cout<<"=================================================================="<<endl;
   cout<<"BH-DFT routine starts here! "<<endl;
   cout<<"Note: "<<endl;
   cout<<"For monometallic clusters: only random xyz moves will be applied "<<endl;
   cout<<"For bimetallic clusters  : 1 atomic swap will be performed after 10 moves "<<endl;
   cout<<"==================================================================="<<endl<<endl;

   i=2;
}

while(i+m < iteraciones)
{
  resto=i%swap_step;
  // Get cluster coordinates from output file
  if(N_Simbolo_2>0) //es bimetalico
  {
    clus_1=extract(Simbolo_1,"geometry.in.next_step");
    clus_2=extract(Simbolo_2,"geometry.in.next_step");
    clus  =clus_1+clus_2;
  }
  else //es monometalico
  {
    clus=extract(Simbolo_1,"geometry.in.next_step");
  }
  // Applies swap or kick
  if(resto==0)
  {
    if(N_Simbolo_1>N_Simbolo_2)
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
      clus.kick(step_width);
    }
    else
    {
      clus.kick_lennard(step_width);
    }
  }
  system("./run.sh");
  contenido=int_pipe("grep 'Have a nice day' output.out | wc -l");
  // Starting randomly if last configuration fails
  while (contenido!=1)
  {
  cout<<" --> SCF failed. Starting again from randomly generated structure! "<<endl;

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
  clus.print_fhi("geometry.tmp");
  system("cat geometry.tmp >> geometry.in ; rm geometry.tmp");
  system("./run.sh");
  contenido=int_pipe("grep 'Have a nice day' output.out | wc -l");
  }
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//                                         SAVE ENERGIES                                          //
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
EnergiaAnterior=Energy;
Energy=float_pipe("grep \" | Total energy of the DFT / Hartree-Fock s.c.f. calculation \" output.out | cut -d \":\" -f 2  | cut -d 'e' -f 1 ");
clus.read_fhi("geometry.in.next_step");
outputfile ="coordinates"+to_string(i)+".xyz";
clus.print_xyz(outputfile);
command.clear(); command="mv output.out output"+to_string(i)+".out";
system(command.c_str());
command.clear(); command="mv geometry.in geometry"+to_string(i)+".in";
system(command.c_str());
command.clear(); command="echo "+to_string(i)+"  "+to_string(Energy)+" >> energies.txt";
system(command.c_str() );
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//                                     Metropolis Monte-Carlo                                     //
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

k_BT = 0.00008617 * Temperature;
if (pow(2.71,(EnergiaAnterior-Energy)/k_BT) > random_number(0,1))
{
  cout<<"--> Basin Hopping MC criteria: Energy accepted! "<<endl;
  cout<<"--> Finished iteration $i"<<endl;
  command.clear(); command="tail -"+to_string(i)+" energies.txt |  sort -nk2 > sorted.txt";
  system(command.c_str());
  i++;
}
else
{
  cout<<"--> Basin Hopping MC criteria: Energy rejected!"<<endl;
  command.clear(); command="mv output"+to_string(i)+".out rejected/outputrejected"+to_string(i)+".out";
  system(command.c_str());
  command.clear(); command="mv geometry"+to_string(i)+".in rejected/geometryrejected"+to_string(i)+".in";
  system(command.c_str());
  m++;
}

}
return 0;
}
