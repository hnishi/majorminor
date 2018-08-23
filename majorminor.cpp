#include"nlib.h"

string majorminor_version = "5.1";
/*
ver2.0 RANGE_COM was implemented
ver3.0 all pairs within 2.0 nm would be written in the output
ver4.0 algorithm was totally replaced --> detection of two minimum contacts between phosphate bacbones
ver4.1 minimum and second minimum --> calculate average of distances of two combination and choose smaller combination
ver4.2 way of residue count was changed to improve judgement of major or minor grooves and bugs were fixed. 
ver4.3 pair was resitricted to i +- 10 bp --> exception was generated 
ver4.4 refactoring
ver5.0 add new mode "base-backbone"
ver5.1 change output format
*/

int flag_debug = 0;

double sq(double i){
   return i*i;
}

int majorminor(  Inp_nishi inp1 ){
   printf("Start majorminor version %s \n", majorminor_version.c_str());

/*
++++++ Input section ++++++
*/   
   string pdbname = inp1.read("INPUTPDB1");
   string outfile1 = inp1.read("OUTPUTFILE1");
   int chain_a_a = atoi(inp1.read("CHAIN_A_a").c_str());
   int chain_a_b = atoi(inp1.read("CHAIN_A_b").c_str());
   int chain_b_a = atoi(inp1.read("CHAIN_B_a").c_str());
   int chain_b_b = atoi(inp1.read("CHAIN_B_b").c_str());
   int chain_c_a = atoi(inp1.read("CHAIN_C_a").c_str());
   int chain_c_b = atoi(inp1.read("CHAIN_C_b").c_str());
   int chain_d_a = atoi(inp1.read("CHAIN_D_a").c_str());
   int chain_d_b = atoi(inp1.read("CHAIN_D_b").c_str());
   double dist_base = atof(inp1.read("DIST_BASE").c_str()); //cutoff distance to define major-backbone mode

   FILE *fout1;
   if((fout1 = fopen(outfile1.c_str(),"w")) == NULL ){
      printf("cannot open output file: %s\n",outfile1.c_str());
      exit(1);
   }

   pdb_nishi* pdb1;
   pdb1 = new pdb_nishi(pdbname.c_str());
   cout<<"TOTAL ATOM = "<< pdb1->total_atom<<endl;
   cout<<"TOTAL RESIDUE = "<< pdb1->total_residue<<endl;

/*
++++++ Search closest point ++++++
Minumum distances betwween phosphorus atoms of all phosphate backbones are calculated
*/
   int pair_residue[10];
   double pair_dist[5] = {999999, 999999, 999999, 999999, 999999};
   int flag_major_backbone; // 0 = major-backbone, 1 = backbone-major

   for(unsigned int i_atom = 0; i_atom < pdb1->total_atom; i_atom++){
      if( ((pdb1->rnum[i_atom] >= chain_a_a && pdb1->rnum[i_atom] <= chain_a_b)
            || (pdb1->rnum[i_atom] >= chain_b_a && pdb1->rnum[i_atom] <= chain_b_b))
            && (pdb1->atmn[i_atom] == "P" 
            || (pdb1->atmn[i_atom] == "O6" && pdb1->resn[i_atom] == "DG")) ){ 
         for(unsigned int j_atom = 0; j_atom < pdb1->total_atom; j_atom++){
            if( ((pdb1->rnum[j_atom] >= chain_c_a && pdb1->rnum[j_atom] <= chain_c_b)
                  || (pdb1->rnum[j_atom] >= chain_d_a && pdb1->rnum[j_atom] <= chain_d_b))
                  && ((pdb1->atmn[i_atom] == "P" && pdb1->atmn[j_atom] == "O6")
                  || (pdb1->atmn[i_atom] == "O6" && pdb1->atmn[j_atom] == "P"))){ 
               double dist = sqrt(sq(pdb1->coox[j_atom] - pdb1->coox[i_atom])
                     +sq(pdb1->cooy[j_atom] - pdb1->cooy[i_atom])
                     +sq(pdb1->cooz[j_atom] - pdb1->cooz[i_atom]));
               if(dist < pair_dist[4]){
                  pair_dist[4] = dist;
                  pair_residue[8] = pdb1->rnum[i_atom];
                  pair_residue[9] = pdb1->rnum[j_atom];
                  if(flag_debug==5)cout<<"DEBUG: pair_dist[4]: "<<pair_dist[4]<<" "<<pair_residue[8]<<" "<<pair_residue[9]<<endl;
                  if(flag_debug==5)cout<<"DEBUG: pdb1->atmn[i_atom], pdb1->atmn[j_atom] = "<<pdb1->atmn[i_atom]<<"-"<<pdb1->atmn[j_atom]<<endl;
                  if( pdb1->atmn[i_atom] == "P" && pdb1->atmn[j_atom] == "O6" ){
                     flag_major_backbone = 1;
                  }
                  else if (pdb1->atmn[i_atom] == "O6" && pdb1->atmn[j_atom] == "P"){
                     flag_major_backbone = 0;
                  }
                  else{return -1;}
               }
            }
         }
      }
      if( pdb1->rnum[i_atom] >= chain_a_a && pdb1->rnum[i_atom] <= chain_a_b 
            && pdb1->atmn[i_atom] == "P" ){ 
         for(unsigned int j_atom = 0; j_atom < pdb1->total_atom; j_atom++){
            if( pdb1->rnum[j_atom] >= chain_c_a && pdb1->rnum[j_atom] <= chain_c_b 
                  && pdb1->atmn[j_atom] == "P" ){ 
               double dist = sqrt(sq(pdb1->coox[j_atom] - pdb1->coox[i_atom])
                     +sq(pdb1->cooy[j_atom] - pdb1->cooy[i_atom])
                     +sq(pdb1->cooz[j_atom] - pdb1->cooz[i_atom]));
               if(dist < pair_dist[0]){
                  pair_dist[0] = dist;
                  pair_residue[0] = pdb1->rnum[i_atom];
                  pair_residue[1] = pdb1->rnum[j_atom];
                  if(flag_debug==1)cout<<"DEBUG: pair_dist[0]: "<<pair_dist[0]<<" "<<pair_residue[0]<<" "<<pair_residue[1]<<endl;
               }
            }
            else if( pdb1->rnum[j_atom] >= chain_d_a && pdb1->rnum[j_atom] <= chain_d_b 
                  && pdb1->atmn[j_atom] == "P" ){ 
               double dist = sqrt(sq(pdb1->coox[j_atom] - pdb1->coox[i_atom])
                     +sq(pdb1->cooy[j_atom] - pdb1->cooy[i_atom])
                     +sq(pdb1->cooz[j_atom] - pdb1->cooz[i_atom]));
               if(dist < pair_dist[1]){
                  pair_dist[1] = dist;
                  pair_residue[2] = pdb1->rnum[i_atom];
                  pair_residue[3] = pdb1->rnum[j_atom];
                  if(flag_debug==1)cout<<"DEBUG: pair_dist[1]: "<<pair_dist[1]<<" "<<pair_residue[2]<<" "<<pair_residue[3]<<endl;
               }
            }
         }
      }
      else if( pdb1->rnum[i_atom] >= chain_b_a && pdb1->rnum[i_atom] <= chain_b_b 
            && pdb1->atmn[i_atom] == "P" ){ 
         for(unsigned int j_atom = 0; j_atom < pdb1->total_atom; j_atom++){
            if( pdb1->rnum[j_atom] >= chain_c_a && pdb1->rnum[j_atom] <= chain_c_b 
                  && pdb1->atmn[j_atom] == "P" ){ 
               double dist = sqrt(sq(pdb1->coox[j_atom] - pdb1->coox[i_atom])
                     +sq(pdb1->cooy[j_atom] - pdb1->cooy[i_atom])
                     +sq(pdb1->cooz[j_atom] - pdb1->cooz[i_atom]));
               if(flag_debug==1)cout<<"DEBUG: dist = "<<dist<<" "<<pdb1->rnum[i_atom]<<" "<<pdb1->rnum[j_atom]<<endl;
               if(flag_debug==1)cout<<"DEBUG: pdb1->coox["<<j_atom<<"], pdb1->coox["<<i_atom<<"]  = "<<pdb1->coox[j_atom]<<", "<< pdb1->coox[i_atom]<<endl;
               if(flag_debug==1)cout<<"DEBUG: pdb1->cooy["<<j_atom<<"], pdb1->cooy["<<i_atom<<"]  = "<<pdb1->cooy[j_atom]<<", "<< pdb1->cooy[i_atom]<<endl;
               if(flag_debug==1)cout<<"DEBUG: pdb1->cooz["<<j_atom<<"], pdb1->cooz["<<i_atom<<"]  = "<<pdb1->cooz[j_atom]<<", "<< pdb1->cooz[i_atom]<<endl;
               if(dist < pair_dist[2]){
                  pair_dist[2] = dist;
                  pair_residue[4] = pdb1->rnum[i_atom];
                  pair_residue[5] = pdb1->rnum[j_atom];
                  if(flag_debug==1)cout<<"DEBUG: pair_dist[2]: "<<pair_dist[2]<<" "<<pair_residue[4]<<" "<<pair_residue[5]<<endl;
               }
            }
            else if( pdb1->rnum[j_atom] >= chain_d_a && pdb1->rnum[j_atom] <= chain_d_b 
                  && pdb1->atmn[j_atom] == "P" ){ 
               double dist = sqrt(sq(pdb1->coox[j_atom] - pdb1->coox[i_atom])
                     +sq(pdb1->cooy[j_atom] - pdb1->cooy[i_atom])
                     +sq(pdb1->cooz[j_atom] - pdb1->cooz[i_atom]));
               if(flag_debug==1)cout<<"DEBUG: zzz = "<<sq(pdb1->coox[j_atom] - pdb1->coox[i_atom])
                     +sq(sq(pdb1->cooy[j_atom] - pdb1->cooy[i_atom]))
                     +sq(sq(pdb1->cooz[j_atom] - pdb1->cooz[i_atom]))<<endl;
               if(flag_debug==1)cout<<"DEBUG: dist = "<<dist<<" "<<pdb1->rnum[i_atom]<<" "<<pdb1->rnum[j_atom]<<endl;
               if(flag_debug==1)cout<<"DEBUG: pdb1->coox["<<j_atom<<"], pdb1->coox["<<i_atom<<"]  = "<<pdb1->coox[j_atom]<<", "<< pdb1->coox[i_atom]<<endl;
               if(flag_debug==1)cout<<"DEBUG: pdb1->cooy["<<j_atom<<"], pdb1->cooy["<<i_atom<<"]  = "<<pdb1->cooy[j_atom]<<", "<< pdb1->cooy[i_atom]<<endl;
               if(flag_debug==1)cout<<"DEBUG: pdb1->cooz["<<j_atom<<"], pdb1->cooz["<<i_atom<<"]  = "<<pdb1->cooz[j_atom]<<", "<< pdb1->cooz[i_atom]<<endl;
               if(dist < pair_dist[3]){
                  pair_dist[3] = dist;
                  pair_residue[6] = pdb1->rnum[i_atom];
                  pair_residue[7] = pdb1->rnum[j_atom];
                  if(flag_debug==1)cout<<"DEBUG: pair_dist[3]: "<<pair_dist[3]<<" "<<pair_residue[6]<<" "<<pair_residue[7]<<endl;
               }
            }
         }
      }
   }



/*
++++++ Which is the mimimum and the second minimum ++++++
*/

   int mode_combination = -1, mode_combination_2 = -1;
   if(pair_dist[0] <= pair_dist[1] && pair_dist[0] <= pair_dist[2] && pair_dist[0] <= pair_dist[3]){
      mode_combination = 0;
      mode_combination_2 = 3;
   }
   else if(pair_dist[1] <= pair_dist[0] && pair_dist[1] <= pair_dist[2] && pair_dist[1] <= pair_dist[3]){
      mode_combination = 1;
      mode_combination_2 = 2;
   }
   else if(pair_dist[2] <= pair_dist[0] && pair_dist[2] <= pair_dist[1] && pair_dist[2] <= pair_dist[3]){
      mode_combination = 2;
      mode_combination_2 = 1;
   }
   else if(pair_dist[3] <= pair_dist[0] && pair_dist[3] <= pair_dist[2] && pair_dist[3] <= pair_dist[1]){
      mode_combination = 3;
      mode_combination_2 = 0;
   }
   else{cerr<<"ERROR: something else"<<endl;}

   if(flag_debug==2)cout<<"DEBUG"<<flag_debug<<": mode_combination = "<<mode_combination<<endl;
   if(flag_debug==2)cout<<"DEBUG"<<flag_debug<<": minimum distance = "<<pair_dist[mode_combination]<<endl;
   if(flag_debug==2)cout<<"DEBUG"<<flag_debug<<": residue pair = "<<pair_residue[mode_combination*2]<<"-"<<pair_residue[mode_combination*2+1]<<endl;
   if(flag_debug==2)cout<<"DEBUG"<<flag_debug<<": mode_combination_2 = "<<mode_combination_2<<endl;
   if(flag_debug==2)cout<<"DEBUG"<<flag_debug<<": second minimum distance = "<<pair_dist[mode_combination_2]<<endl;
   if(flag_debug==2)cout<<"DEBUG"<<flag_debug<<": residue pair = "<<pair_residue[mode_combination_2*2]<<"-"<<pair_residue[mode_combination_2*2+1]<<endl;


/*
++++++ Calculate average of distances of two combinations ++++++
ver 4.1
*/
/* 
   int mode_combination = -1, mode_combination_2 = -1;
   double ave[2];
   ave[0] = (pair_dist[0] + pair_dist[3])/2;
   ave[1] = (pair_dist[1] + pair_dist[2])/2;
   if(ave[0] < ave[1]){
      mode_combination = 0;
      mode_combination_2 = 3;
   }
   else{
      mode_combination = 1;
      mode_combination_2 = 2;
   }
*/   

/*
++++++ Judgement of major or minor grooves ++++++
*/
   string mode_majorminor_1, mode_majorminor_2; 
   int tmp_resi1 = pair_residue[mode_combination*2] - chain_a_a + 1;
   //int tmp_resi2 = pair_residue[mode_combination_2*2] - chain_b_a + 1; # before ver4.2
   int tmp_resi2 = chain_b_b - pair_residue[mode_combination_2*2] + 1;
   if( sqrt(sq(tmp_resi1 - tmp_resi2)) > 10 ){
      mode_majorminor_1 = "exception";
      cerr<<"ERROR: "<<pdbname<<endl;
      cerr<<"ERROR: Major or minor grooves cannot be defined."<<endl;
      cerr<<"ERROR: Difference between two nucleotides is too large (> 10). "<<endl;
      cerr<<"ERROR: tmp_resi1 = "<<tmp_resi1<<endl;
      cerr<<"ERROR: tmp_resi2 = "<<tmp_resi2<<endl;
   }
   else if( tmp_resi1 > tmp_resi2 ){
      mode_majorminor_1 = "minor";
   }
   else if( tmp_resi1 < tmp_resi2 ){
      mode_majorminor_1 = "major";
   }
   else{
      cerr<<"ERROR: "<<pdbname<<endl;
      cerr<<"ERROR: Major or minor grooves cannot be defined."<<endl;
      cerr<<"ERROR: tmp_resi1 = "<<tmp_resi1<<endl;
      cerr<<"ERROR: tmp_resi2 = "<<tmp_resi2<<endl;
   }

   int tmp_resi3 = pair_residue[mode_combination*2+1] - chain_c_a + 1;
   //int tmp_resi4 = pair_residue[mode_combination_2*2+1] - chain_d_a + 1; # before ver4.2
   int tmp_resi4 = chain_d_b - pair_residue[mode_combination_2*2+1] + 1;
   if( sqrt(sq(tmp_resi3 - tmp_resi4)) > 10 ){
      mode_majorminor_2 = "exception";
      cerr<<"ERROR: "<<pdbname<<endl;
      cerr<<"ERROR: Major or minor grooves cannot be defined."<<endl;
      cerr<<"ERROR: Difference between two nucleotides is too large (> 10). "<<endl;
      cerr<<"ERROR: tmp_resi3 = "<<tmp_resi3<<endl;
      cerr<<"ERROR: tmp_resi4 = "<<tmp_resi4<<endl;
   }
   else if( tmp_resi3 > tmp_resi4 ){
      mode_majorminor_2 = "minor";
   }
   else if( tmp_resi3 < tmp_resi4 ){
      mode_majorminor_2 = "major";
   }
   else{
      cerr<<"ERROR: "<<pdbname<<endl;
      cerr<<"ERROR: Major or minor grooves cannot be defined."<<endl;
      cerr<<"ERROR: tmp_resi3 = "<<tmp_resi3<<endl;
      cerr<<"ERROR: tmp_resi4 = "<<tmp_resi4<<endl;
   }
   if( pair_dist[4] < dist_base ){
      if( flag_major_backbone == 0 ){ // 0 = major-backbone, 1 = backbone-major
         mode_majorminor_1 = "base";
         mode_majorminor_2 = "backbone";
      }
      else if( flag_major_backbone == 1 ){
         mode_majorminor_1 = "backbone";
         mode_majorminor_2 = "base";
      }
   }
   cout<<"Mode = "<<mode_majorminor_1<<"-"<<mode_majorminor_2<<endl;

/*
++++++ Output file ++++++
*/
/*   if( pair_dist[4] < dist_base ){
      fprintf(fout1,"%s-%s %d-%d %.3f %s \n", 
         mode_majorminor_1.c_str(), mode_majorminor_2.c_str(), 
         pair_residue[8], pair_residue[9], pair_dist[4], pdbname.c_str());
   }
   else{
      fprintf(fout1,"%s-%s %d-%d %.3f %d-%d %.3f %s \n", 
         mode_majorminor_1.c_str(), mode_majorminor_2.c_str(), 
         pair_residue[mode_combination*2], pair_residue[mode_combination*2+1], pair_dist[mode_combination],
         pair_residue[mode_combination_2*2], pair_residue[mode_combination_2*2+1], pair_dist[mode_combination_2],
         pdbname.c_str());
   }
*/
   fprintf(fout1,"%s-%s %d-%d %.3f %d-%d %.3f %d-%d %.3f %s \n", 
      mode_majorminor_1.c_str(), mode_majorminor_2.c_str(), 
      pair_residue[mode_combination*2], pair_residue[mode_combination*2+1], pair_dist[mode_combination],
      pair_residue[mode_combination_2*2], pair_residue[mode_combination_2*2+1], pair_dist[mode_combination_2],
      pair_residue[8], pair_residue[9], pair_dist[4], pdbname.c_str());

   fclose(fout1);
   return 0;
}

