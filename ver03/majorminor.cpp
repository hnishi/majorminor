#include"nlib.h"

/*
VERSION majorminor 3.0

ver2.0 RANGE_COM was implemented
ver3.0 all pairs within 2.0 nm would be written in the output
*/

double sq(double i){
   return i*i;
}

int majorminor(  Inp_nishi inp1 ){
   puts("Start majorminor");
   
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
   int range_com = atoi(inp1.read("RANGE_COM").c_str());
   double cutoff_dist = atof(inp1.read("CUTOFF_DIST").c_str());
/*
   cout<<"DEBUG: "<<chain_a_a<<endl;
   cout<<"DEBUG: "<<chain_a_b<<endl;
   cout<<"DEBUG: "<<chain_b_a<<endl;
   cout<<"DEBUG: "<<chain_b_b<<endl;
*/

   FILE *fout1;
   if((fout1 = fopen(outfile1.c_str(),"w")) == NULL ){
      printf("cannot open output file: %s\n",outfile1.c_str());
      exit(1);
   }

   pdb_nishi* pdb1;
   pdb1 = new pdb_nishi(pdbname.c_str());
   cout<<"TOTAL ATOM = "<< pdb1->total_atom<<endl;
   cout<<"TOTAL RESIDUE = "<< pdb1->total_residue<<endl;


   vector<double> vec_com_major_dsDNA1; // 3N dimension; i.e. (x,y,z)*N
   vector<double> vec_com_major_dsDNA2;
   int ii=0;

   ii = chain_b_b - range_com;
   for(int i=chain_a_a + range_com;i<=chain_a_b - range_com;i++){
      double comx=0,comy=0,comz=0;
      int n_factor=0;
      //cout<<i<<endl;
      //cout<<ii<<endl;
      for(unsigned int j=0;j<pdb1->total_atom;j++){
         //cout<<"DEBUG: rnum "<<pdb1->rnum[j]<<endl;
         if( ( i >= chain_a_a + range_com && i <= chain_a_b - range_com )
            && ( ii >= chain_b_a + range_com && ii <= chain_b_b - range_com )
            && ( ( pdb1->rnum[j] >= i - range_com && pdb1->rnum[j] <= i + range_com ) 
               || ( pdb1->rnum[j] >= ii - range_com && pdb1->rnum[j] <= ii + range_com ) )
            && (pdb1->atmn[j] == "N1" || pdb1->atmn[j] == "C2"
               || pdb1->atmn[j] == "N3"
               || pdb1->atmn[j] == "C4"
               || pdb1->atmn[j] == "C5"
               || pdb1->atmn[j] == "C6"
               || pdb1->atmn[j] == "N7"
               || pdb1->atmn[j] == "C8"
               || pdb1->atmn[j] == "N9" ) ){
            //cout<<"hit!!! "<<pdb1->rnum[j]<<pdb1->atmn[j]<<" "<<pdb1->cooz[j]<<endl;
            comx=comx+pdb1->coox[j];
            comy=comy+pdb1->cooy[j];
            comz=comz+pdb1->cooz[j];
            n_factor++;
         }
      }
      comx=comx/double(n_factor);
      comy=comy/double(n_factor);
      comz=comz/double(n_factor);
      vec_com_major_dsDNA1.push_back(comx);
      vec_com_major_dsDNA1.push_back(comy);
      vec_com_major_dsDNA1.push_back(comz);
      cout<<"DEBUG: residue= "<<i<<" "<<ii<<endl;
      cout<<"DEBUG: n_factor= "<<n_factor<<endl;
      cout<<"DEBUG: comx= "<<comx<<endl;
      cout<<"DEBUG: comy= "<<comy<<endl;
      cout<<"DEBUG: comz= "<<comz<<endl;
      ii--;
   }
   cout<<"DEBUG: "<<vec_com_major_dsDNA1.size()<<endl;

   ii = chain_d_b - range_com;
   for(int i=chain_c_a + range_com ;i<=chain_c_b - range_com;i++){
      double comx=0,comy=0,comz=0;
      int n_factor=0;
      //cout<<i<<endl;
      //cout<<ii<<endl;
      for(unsigned int j=0;j<pdb1->total_atom;j++){
         //cout<<"DEBUG: rnum "<<pdb1->rnum[j]<<endl;
         if( ( i >= chain_c_a + range_com && i <= chain_c_b - range_com )
            && ( ii >= chain_d_a + range_com && ii <= chain_d_b - range_com )
            && ( ( pdb1->rnum[j] >= i - range_com && pdb1->rnum[j] <= i + range_com ) 
               || ( pdb1->rnum[j] >= ii - range_com && pdb1->rnum[j] <= ii + range_com ) )
            && (pdb1->atmn[j] == "N1" || pdb1->atmn[j] == "C2"
               || pdb1->atmn[j] == "N3"
               || pdb1->atmn[j] == "C4"
               || pdb1->atmn[j] == "C5"
               || pdb1->atmn[j] == "C6"
               || pdb1->atmn[j] == "N7"
               || pdb1->atmn[j] == "C8"
               || pdb1->atmn[j] == "N9" ) ){
            //cout<<"hit!!! "<<pdb1->rnum[j]<<pdb1->atmn[j]<<" "<<pdb1->cooz[j]<<endl;
            comx=comx+pdb1->coox[j];
            comy=comy+pdb1->cooy[j];
            comz=comz+pdb1->cooz[j];
            n_factor++;
         }
      }
      comx=comx/double(n_factor);
      comy=comy/double(n_factor);
      comz=comz/double(n_factor);
      vec_com_major_dsDNA2.push_back(comx);
      vec_com_major_dsDNA2.push_back(comy);
      vec_com_major_dsDNA2.push_back(comz);
      cout<<"DEBUG: n_factor= "<<n_factor<<endl;
      cout<<"DEBUG: comx= "<<comx<<endl;
      cout<<"DEBUG: comy= "<<comy<<endl;
      cout<<"DEBUG: comz= "<<comz<<endl;
      ii--;
   }
   cout<<"DEBUG: "<<vec_com_major_dsDNA2.size()<<endl;

   vector<double> vec_com_minor_dsDNA1;
   vector<double> vec_com_minor_dsDNA2;

   ii = chain_b_b;
   for(int i=chain_a_a;i<=chain_a_b;i++){
      double comx=0,comy=0,comz=0;
      int n_factor=0;
      //cout<<i<<endl;
      //cout<<ii<<endl;
      for(unsigned int j=0;j<pdb1->total_atom;j++){
         //cout<<"DEBUG: rnum "<<pdb1->rnum[j]<<endl;
         if( (pdb1->rnum[j] == i || pdb1->rnum[j] == ii )
            && (pdb1->atmn[j] == "C1'" || pdb1->atmn[j] == "C2'"
            || pdb1->atmn[j] == "C3'"
            || pdb1->atmn[j] == "C4'"
            || pdb1->atmn[j] == "C5'"
            || pdb1->atmn[j] == "O4'" ) ){
            //cout<<"hit!!! "<<pdb1->rnum[j]<<pdb1->atmn[j]<<" "<<pdb1->cooz[j]<<endl;
            comx=comx+pdb1->coox[j];
            comy=comy+pdb1->cooy[j];
            comz=comz+pdb1->cooz[j];
            n_factor++;
         }
      }
      comx=comx/double(n_factor);
      comy=comy/double(n_factor);
      comz=comz/double(n_factor);
      vec_com_minor_dsDNA1.push_back(comx);
      vec_com_minor_dsDNA1.push_back(comy);
      vec_com_minor_dsDNA1.push_back(comz);
      cout<<"DEBUG: n_factor= "<<n_factor<<endl;
      cout<<"DEBUG: comx= "<<comx<<endl;
      cout<<"DEBUG: comy= "<<comy<<endl;
      cout<<"DEBUG: comz= "<<comz<<endl;
      ii--;
   }
   cout<<"DEBUG: "<<vec_com_minor_dsDNA1.size()<<endl;

   ii = chain_d_b;
   for(int i=chain_c_a;i<=chain_c_b;i++){
      double comx=0,comy=0,comz=0;
      int n_factor=0;
      //cout<<i<<endl;
      //cout<<ii<<endl;
      for(unsigned int j=0;j<pdb1->total_atom;j++){
         //cout<<"DEBUG: rnum "<<pdb1->rnum[j]<<endl;
         if( (pdb1->rnum[j] == i || pdb1->rnum[j] == ii )
            && (pdb1->atmn[j] == "C1'" || pdb1->atmn[j] == "C2'"
            || pdb1->atmn[j] == "C3'"
            || pdb1->atmn[j] == "C4'"
            || pdb1->atmn[j] == "C5'"
            || pdb1->atmn[j] == "O4'" ) ){
            //cout<<"hit!!! "<<pdb1->rnum[j]<<pdb1->atmn[j]<<" "<<pdb1->cooz[j]<<endl;
            comx=comx+pdb1->coox[j];
            comy=comy+pdb1->cooy[j];
            comz=comz+pdb1->cooz[j];
            n_factor++;
         }
      }
      comx=comx/double(n_factor);
      comy=comy/double(n_factor);
      comz=comz/double(n_factor);
      vec_com_minor_dsDNA2.push_back(comx);
      vec_com_minor_dsDNA2.push_back(comy);
      vec_com_minor_dsDNA2.push_back(comz);
      cout<<"DEBUG: n_factor= "<<n_factor<<endl;
      cout<<"DEBUG: comx= "<<comx<<endl;
      cout<<"DEBUG: comy= "<<comy<<endl;
      cout<<"DEBUG: comz= "<<comz<<endl;
      ii--;
   }
   cout<<"DEBUG: "<<vec_com_minor_dsDNA2.size()<<endl;

   //for(unsigned int i=0;i<pdb1->total_atom;i++){
   //   printf("(i,x,y,x)=(%d,%.3lf,%.3lf,%.3lf) \n",pdb1->anum[i],pdb1->coox[i],pdb1->cooy[i],pdb1->cooz[i]);
   //}


   string mode = "none";
   int bp_dsDNA1 = -1, bp_dsDNA2 = -1;
   double min_dist = 999999,dist;

   for(unsigned int i=0;i<vec_com_major_dsDNA1.size();i=i+3 ){
      //cout<<"DEBUG: "<<i<<endl;
      for(unsigned int j=0;j<vec_com_major_dsDNA2.size();j=j+3 ){
         dist = sq(vec_com_major_dsDNA1[i]-vec_com_major_dsDNA2[j])
                + sq(vec_com_major_dsDNA1[i+1]-vec_com_major_dsDNA2[j+1])
                + sq(vec_com_major_dsDNA1[i+2]-vec_com_major_dsDNA2[j+2]);
         //cout<<"DEBUG: dist = "<<sqrt(dist)<<endl;
         bp_dsDNA1 = i/3+1+range_com+chain_a_a;
         bp_dsDNA2 = j/3+1+range_com+chain_c_a;
         mode = "major-major" ; 
         if( dist < min_dist ){
            min_dist = dist;
            cout<<"DEBUG: min_dist = "<<sqrt(min_dist)<<endl;
            cout<<"DEBUG: mode = "<<mode<<endl;
            cout<<"DEBUG: bp_DNA1 = "<<bp_dsDNA1<<endl;
            cout<<"DEBUG: bp_DNA2 = "<<bp_dsDNA2<<endl;
         }
         if( sqrt(dist) < cutoff_dist ){
            fprintf(fout1,"%s %d-%d %.3f %s \n", mode.c_str(),bp_dsDNA1,bp_dsDNA2,sqrt(dist),pdbname.c_str());
         } 
      }
      for(unsigned int j=0;j<vec_com_minor_dsDNA2.size();j=j+3 ){
         dist = sq(vec_com_major_dsDNA1[i]-vec_com_minor_dsDNA2[j])
                + sq(vec_com_major_dsDNA1[i+1]-vec_com_minor_dsDNA2[j+1])
                + sq(vec_com_major_dsDNA1[i+2]-vec_com_minor_dsDNA2[j+2]);
         //cout<<"DEBUG: dist = "<<sqrt(dist)<<endl;
         bp_dsDNA1 = i/3+1+range_com+chain_a_a;
         bp_dsDNA2 = j/3+1+range_com+chain_c_a;
         mode = "major-minor" ; 
         if( dist < min_dist ){
            min_dist = dist;
            //cout<<"DEBUG: min_dist = "<<sqrt(min_dist)<<endl;
         }
         if( sqrt(dist) < cutoff_dist ){
            fprintf(fout1,"%s %d-%d %.3f %s \n", mode.c_str(),bp_dsDNA1,bp_dsDNA2,sqrt(dist),pdbname.c_str());
         } 
      }
   }
   for(unsigned int i=0;i<vec_com_minor_dsDNA1.size();i=i+3 ){
      //cout<<"DEBUG: "<<i<<endl;
      for(unsigned int j=0;j<vec_com_major_dsDNA2.size();j=j+3 ){
         dist = sq(vec_com_minor_dsDNA1[i]-vec_com_major_dsDNA2[j])
                + sq(vec_com_minor_dsDNA1[i+1]-vec_com_major_dsDNA2[j+1])
                + sq(vec_com_minor_dsDNA1[i+2]-vec_com_major_dsDNA2[j+2]);
         //cout<<"DEBUG: dist = "<<sqrt(dist)<<endl;
         bp_dsDNA1 = i/3+1+range_com+chain_a_a;
         bp_dsDNA2 = j/3+1+range_com+chain_c_a;
         mode = "minor-major" ; 
         if( dist < min_dist ){
            min_dist = dist;
            //cout<<"DEBUG: min_dist = "<<sqrt(min_dist)<<endl;
         }
         if( sqrt(dist) < cutoff_dist ){
            fprintf(fout1,"%s %d-%d %.3f %s \n", mode.c_str(),bp_dsDNA1,bp_dsDNA2,sqrt(dist),pdbname.c_str());
         } 
      }
      for(unsigned int j=0;j<vec_com_minor_dsDNA2.size();j=j+3 ){
         dist = sq(vec_com_minor_dsDNA1[i]-vec_com_minor_dsDNA2[j])
                + sq(vec_com_minor_dsDNA1[i+1]-vec_com_minor_dsDNA2[j+1])
                + sq(vec_com_minor_dsDNA1[i+2]-vec_com_minor_dsDNA2[j+2]);
         //cout<<"DEBUG: dist = "<<sqrt(dist)<<endl;
         bp_dsDNA1 = i/3+1+range_com+chain_a_a;
         bp_dsDNA2 = j/3+1+range_com+chain_c_a;
         mode = "minor-minor" ; 
         if( dist < min_dist ){
            min_dist = dist;
            cout<<"DEBUG: min_dist = "<<sqrt(min_dist)<<endl;
         }
         if( sqrt(dist) < cutoff_dist ){
            fprintf(fout1,"%s %d-%d %.3f %s \n", mode.c_str(),bp_dsDNA1,bp_dsDNA2,sqrt(dist),pdbname.c_str());
         } 
      }
   }
   //fprintf(fout1,"mode pair distance file \n");
   //fprintf(fout1,"%s %d-%d %.3f %s \n", mode.c_str(),bp_dsDNA1,bp_dsDNA2,sqrt(min_dist),pdbname.c_str());


   fclose(fout1);
   return 0;
}
