#include<atomicpp.h>

string Simbolo_1, Simbolo_2;
int N_Simbolo_1, N_Simbolo_2, n;
string initialization_file;
int continue_alg,  Ncore, randomness, kick, iteraciones,swap_step, contenido, m;
float step_width, Temperature ;
string file_name, command;

int main(int argc, char *argv[])
{
	
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//                                    Gets data from input.bh                                     //
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/


Simbolo_1=string_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 2 | cut -d ':' -f 1 ");
Simbolo_2=string_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 3 | cut -d ':' -f 1 ");
N_Simbolo_1=int_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 2 | cut -d ':' -f 2 | cut -d ']' -f 1 ");
N_Simbolo_2=int_pipe("grep 'cluster_ntyp' input.bh | cut -d '[' -f 3 | cut -d ':' -f 2 | cut -d ']' -f 1 ");
continue_alg=int_pipe("grep 'continue' input.bh | cut -d ' ' -f 3 >> variables");
initialization_file=string_pipe("grep 'initialization_file' input.bh | awk '{print $3}' ");
randomness=int_pipe("grep 'randomness' input.bh | awk '{print $3}' ");
kick=int_pipe("grep 'kick_type' input.bh | awk '{print $3}' ");
file_name=string_pipe("grep 'file_name' input.bh | awk '{print $3}' ");
step_width=float_pipe("grep 'step_width' input.bh | awk '{print $3}' "); //o double?
Temperature=float_pipe("grep 'temperature_K' input.bh | awk '{ print $3 }' "); //o double?
Ncore=int_pipe("grep 'Ncore' input.bh | head -1 | awk '{print $3}' ");
iteraciones=system("grep 'iterations' input.bh | awk '{ print $3 }' ");
//m=int_pipe("......");
swap_step=system("grep 'swap_step' input.bh | awk '{ print $3 }' ");

cout<<Simbolo_1<<" "<<Simbolo_2<<" "<<N_Simbolo_1<<" "<<N_Simbolo_2<<" "<<n<<" "<<continue_alg<<" "<<initialization_file<<" "<<randomness<<" "<<kick<<" "<<file_name<<" "<<step_width<<" "<<Temperature<<"  "<<iteraciones<<" "<<swap_step<<endl;
//################################# CREATES Work Directory #####################################
command="if [ -d $file_name ] ; then mv $file_name other_$file_name ; fi ; cp -r input $file_name";
system(comand);

//##############################################################################################
//#                                         BEGIN ALGORITHM                                    #
//##############################################################################################

system("cd $file_name ; mkdir rejected");
contenido=0;
count=1;
while(contenido!=1)
{

  if(initialization_file!=false && count==1)
  {
    //Inicializa archivo geometry.in
    system("cp initalization_file geometry.in")
    count++;
  }
  else
  {
    //Empieza desde 0
    if(type == "bimetalic")  //Esto quiere decir que es bimetalico
    {
      if(randomness==1)  // fully random
      {
        clus.srand_generator
      }
      else //pseudoaleatorio
      {
        clus.rand_generator
      }
    }
    else //es monometalico
    {
    	if(randomness==1)  // fully random
      {
        clus.srand_generator
      }
      else //pseudoaleatorio
      {
        clus.rand_generator
      }
    }
  }
    
  system("./run.sh");
  system("grep 'Have a nice day.' output.out | wc -l > contenido.txt"); //Checar output nombre, poner if=1 then
  ifstream cont("contenido.txt");
  cont>>contenido;
  cont.close();
  system("rm contenido.txt");

}

system("grep \" | Total energy of the DFT / Hartree-Fock s.c.f. calculation \" output.out | cut -d \":\" -f 2 > energy.txt")
ifstream en("energy.txt");
en>>energy;
en.close();
system("rm energy.txt");

clus.read_fhi("geometry.in.next_step");
clus.print_xyz("coordinates1.xyz");
system("mv output.out output1.out");
system("mv geometry.in geometry1.in");
system("echo Step Energy >> energies.txt");
system("echo 1     "+to_string(energy)+" >> energies.txt".c_str() );

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
    clus.swap(swap_step);
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
      clus.kicklenard(step_width);
    }
  }
  system("./run.sh");
  system("grep 'Have a nice day. ' output.out | wc -l > contenido.txt");
  ifstream cont("contenido.txt");
  cont>>contenido;
  cont.close();
  system("rm contenido.txt");
  

  while (contenido!=1)
  {
  cout<<" --> SCF failed. Starting again from randomly generated structure! "<<endl;

  if(n>3)  //Esto quiere decir que es bimetalico
  {
    	if(randomness==1)  // fully random
      {
        clus.srand_generator
      }
      else //pseudoaleatorio
      {
        clus.rand_generator
      }
  }
  else //es monometalico
  {
      	if(randomness==1)  // fully random
      {
        clus.srand_generator
      }
      else //pseudoaleatorio
      {
        clus.rand_generator
      }
  }

  system("./run.sh");
  system("grep 'Have a nice day. ' output.out | wc -l > contenido.txt");
  ifstream cont("contenido.txt");
  cont>>contenido;
  cont.close();
  system("rm contenido.txt");
  }

  ################################################################################################
  #                                     Save  configuration                                      #
  ################################################################################################


EnergiaAnterior=energy;
system("grep \" | Total energy of the DFT / Hartree-Fock s.c.f. calculation \" output.out | cut -d \":\" -f 2  | cut -d 'e' -f 1 > energy.txt");
ifstream en("energy.txt");
en>>energy;
en.close();
system("rm energy.txt");

clus.read_fhi("geometry.in.next_step");
clus.print_xyz("cordinates"+to_string(i)+".xyz".c_str(), energy );
system("mv output.out output"+to_string(i)+".out".c_str());
system("mv geometry.in geometry"+to_string(i)+".in".c_str());
system("echo "+to_string(i)+"  "+to_string(energy)+" >> energies.txt".c_str() );


  ################################################################################################
  #                                  Metropolis Monte-Carlo                                      #
  ################################################################################################


kBT = 0.00008617 * temperature_K;
if (pow(e,(EnergiaAnterior-energy)/kBT) > random_number(0,1))
{
	cout<<"--> Basin Hopping MC criteria: Energy accepted! "<<endl;
    cout<<"--> Finished iteration $i"<<endl;
    system("tail "-$i "energies.txt |  sort -nk2 > sorted.txt");
    i++;
}
else
{
	cout<<echo "--> Basin Hopping MC criteria: Energy rejected!"<<endl;
    system("mv output"$i".out rejected/outputrejected"$i".out");
    system("mv geometry"$i".in rejected/geometryrejected"$i".in");
    m++;
}
}
return 0;
}
