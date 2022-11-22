// [[Rcpp::depends(qs)]]
#include <Rcpp.h>
#include <qs.h>
using namespace Rcpp;
// [[Rcpp::plugins(unwindProtect)]]
//unwind protect (requires R>3.5.x), plugin speeds up use of R functions inside Rcpp



//functions that takes GRbaseCoverage and output GRcoverageInbins output
//taken from GRcoverageInbins function from compEpiTools
// and re-implemented in Rcpp


//no chunks. should be efficient enough
// [[Rcpp::export]]
List encryptcov (List coverage){
  //coverage: is a list of integerVectors from coverage
  int n = coverage.size();
  //define output arrays
  Rcpp::IntegerVector keys(n);
  Rcpp::List enc_cov(n);
  Rcpp::List final_list(2);

  //here define all the keys
  //1- 2-byte encoding
  //2- 1-byte encoding
  //3- nibble encoding, even
  //4- nibble encoding, odd
  
  //loop through the list. Keep condition
  for (int i=0; i<n; i++ ){
    Rcpp::IntegerVector currentrange(coverage[i]);
    int arraylength=currentrange.size();
    
    int maxElement = *std::max_element(currentrange.begin(), currentrange.end());

    if(maxElement<=15){
      //int compressedval;
      int counter=0;
      // 1/2 byte (nibble) compression. Separate even and odd
      if(arraylength%2 ==0){

        //even
        Rcpp::RawVector currentencrypted(arraylength/2);
        keys[i]=3;
        for (int j=0; j<arraylength; j+=2){
          currentencrypted[counter]=currentrange[j]*16+currentrange[j+1];
          counter++;
        }
        enc_cov[i]=currentencrypted;

      }else{

        //odd
        Rcpp::RawVector currentencrypted(arraylength/2+1);
        keys[i]=4;
        for (int j=0; j<arraylength-1; j+=2){
          currentencrypted[counter]=currentrange[j]*16+currentrange[j+1];
          counter++;
        }
        currentencrypted[counter]=currentrange[arraylength-1];
        enc_cov[i]=currentencrypted;

      }

    }else if (maxElement<=255){

      //1-byte compression; just change objcet type without loop!
      keys[i]=2;
      Rcpp::RawVector currentencrypted(arraylength);
      currentencrypted=currentrange;
      enc_cov[i]=currentencrypted;

    }else{

      //number is >255, 2-byte compression
      keys[i]=1;
      Rcpp::RawVector currentencrypted(arraylength*2);
      int counter=0;
      for (int j=0; j<arraylength; j++){
        currentencrypted[counter]=currentrange[j]/256;
        currentencrypted[counter+1]=currentrange[j]%256;
        counter+=2;
      }
      enc_cov[i]=currentencrypted;

    }

  }


  final_list[0]=enc_cov;
  final_list[1]=keys;
  
  return(final_list);
}





//cut the transcripts ranges in specific indexes and sum the coverage inside the internal part (remaining range),
//the input is always a list of kind "GRbaseCoverageOutput", from GRbaseCoverage2 function
//StartingPositions and EndingPositions are the array of positions for each element of the GRbaseCoverageOutput list
//it returns an integer vector, that is the sums of the elements of the list within the cut part (be careful to the 0-index of C!)

// [[Rcpp::export]]
NumericVector cutAndSumTranscriptsCPP (List GRbaseCoverageOutput,IntegerVector StartingPositions,IntegerVector EndingPositions, IntegerVector keys){  
    //Rcpp::List xlist(GRbaseCoverageOutput); 
    //Rcpp::IntegerVector starts(StartingPositions);
    //Rcpp::IntegerVector ends(EndingPositions);

    int n = GRbaseCoverageOutput.size(); 
    //std::vector<int> lengths(n);
    int res;
    //define final result
    //std::vector<double> results(n);
    Rcpp::NumericVector results(n);

    Function decompress ("lz4_decompress_raw");
    for(int i=0; i<n; i++) {   


      Rcpp::RawVector currentrange=decompress(GRbaseCoverageOutput[i]);
      //Rcpp::RawVector currentrange(GRbaseCoverageOutput[i]);  
      int currentkeychoice=keys[i];

      //HERE DIFFERENT KIND OF COMPRESSIONS for each range
      //1-byte compression (=FALSE)
      if (currentkeychoice==2){
        res=0;
        // for each interval, sum
        for(int k=StartingPositions[i]; k< EndingPositions[i]; k++){
          //printf("cumulative: %f\\n",res);
          //sum the cut element
          res+=currentrange[k]; 
        }
        //populate the final array with the cut sum
        results[i]=res;        
      }


      //nibble compression even/odd
      if (currentkeychoice==3 | currentkeychoice==4){
        //for each interval, decrypt and sum
        int currenthalf;
        int counterblock=StartingPositions[i]/2; //+1??
        //if no reminder, start from the most signif. nibble of the next byte, otherwise from the least signif.
        if(StartingPositions[i]%2==1){
          currenthalf=1;
        }else{
          currenthalf=0;
        }

        unsigned char num_byte;
        int res=0;
        // for each bin, sum elements in each bin. The remaining is considered 
        for(int k=StartingPositions[i]; k< EndingPositions[i]; k++){
          if(currenthalf==0  ){
            num_byte= (currentrange[counterblock] & 0xF0) >>4 ;
          //else, take only last 4 bits (second nibble)
          }else{
            num_byte= currentrange[counterblock] & 0x0F ;     
          }
          res+=num_byte;


          currenthalf++;
          if (currenthalf==2){
            //reset, this byte expired
            counterblock++;
            currenthalf=0;
          }
          //if 4 (odd) and nbin-1, add last byte of the currentrange (otherwise only most significant nibble is added (=0 by definition))

          // if(currentkeychoice==4 & counterblock==m-1 ){
          //   num_byte=currentrange[m-1] & 0x0F;
          //   res+=num_byte;
          // }  

        }
        results[i]=res;
      }



      //2-byte compression (=TRUE)
      if (currentkeychoice==1){
        //loop through the elements of the current range, multiplied by 2
        int temp=0;
        //for each chunk (little array), sum everyhing 2 by 2 (2-byte compression, big endian)
        for(int k=StartingPositions[i]; k<EndingPositions[i]; k+=2){     
          temp+=currentrange[k]*256+currentrange[k+1];
        }   
        results[i]=temp;                 
      }


      //lengths[i]=m;

    }
    return (results); 
}


//thanks to SymbolixAU (which_maxCpp3 function) help in https://stackoverflow.com/questions/59055902/find-index-of-all-max-min-values-in-vector-in-rcpp
//It works, but it takes the FIRST max value and not a random one for each range
// [[Rcpp::export]]
IntegerVector findPosSummitCPP(List GRbaseCoverageOutput, IntegerVector keys){
  int n = GRbaseCoverageOutput.size();
  Rcpp::IntegerVector finalpositions(n);
  int pos_max;
  Function decompress ("lz4_decompress_raw");
  for(int i=0; i<n; i++) {
    Rcpp::RawVector currentrange=decompress(GRbaseCoverageOutput[i]);
    //Rcpp::RawVector currentrange(GRbaseCoverageOutput[i]);
    int currentkeychoice=keys[i];
    int m=currentrange.size();

    //pos_max already 1-based for R
    pos_max=1;

    //HERE DIFFERENT KIND OF COMPRESSIONS for each range

    //1-byte compression 
    if (currentkeychoice==2){
      int current_max=currentrange[0];
      //start from 1 (0 element is the starting one)
      for (int j=1; j<m; j++){
        int currentval=currentrange[j];
        //new max. update everything
        if(currentval>current_max){
          pos_max=j+1;
          current_max=currentval;
        //same max. need to keep and then random choice ?
        }else if (currentval==current_max){
          //
        }
      }  
    }


    //nibble compression even/odd
    if (currentkeychoice==3 | currentkeychoice==4){
      int current_max=0;
      int iterations;
      //for each interval, decrypt and sum
      if (currentkeychoice==3){
        iterations=m*2;
      }else{
        iterations=m*2-1;
      }
      int currenthalf=0;
      int counterblock=0;
      unsigned char num_byte;
      int currentval;
      // for each bin, sum elements in each bin. The remaining is considered 
      for(int k=0; k< iterations; k++){
        if(currenthalf==0  ){
          num_byte= (currentrange[counterblock] & 0xF0) >>4 ;
        //else, take only last 4 bits (second nibble)
        }else{
          num_byte= currentrange[counterblock] & 0x0F ;     
        }
        currentval=num_byte;
        
        if (currentval>current_max){
          pos_max=k+1;
          current_max=currentval;
        }

        currenthalf++;
        if (currenthalf==2){
          //reset, this byte expired
          counterblock++;
          currenthalf=0;
        }
        //if 4 (odd) and nbin-1, add last byte of the currentrange (otherwise only most significant nibble is added (=0 by definition))
        //AND only if Nbins can be divided into chunks?? If remaining, do not do it
        if(currentkeychoice==4 & k==iterations-1 ){
          num_byte=currentrange[m-1] & 0x0F;
          currentval=num_byte;
          if(currentval>current_max){
            pos_max=k+2;
          }
        }          
      }

    }
    
    
    //2-byte compression 
    if(currentkeychoice==1){
      int counter=1;
      int temp;
      int current_max=currentrange[0]*256+currentrange[1];
      for(int j=0; j<m; j+=2){     
        temp=currentrange[j]*256+currentrange[j+1];
        //new max. update everything
        if(temp>current_max){
          pos_max=counter;
          current_max=temp;
        //same max. need to keep and then random choice ?
        }else if (temp==current_max){
          //
        }        
        counter++;
      } 
    }


    finalpositions[i]=pos_max;

  }
  return(finalpositions);
}







//extremely efficient implementation of zero check of coverages in CPP. No more RAM peaks
//much faster

// [[Rcpp::export]]
LogicalVector verifyzerocov(List covresult){  
     //Rcpp::List xlist(covresult); 

     int n = covresult.size(); 
     //bool allzeros=true;
     //loop through all ranges in list
     int tempsum=0;
     Function decompress ("lz4_decompress_raw");
     Rcpp::LogicalVector ALLzeros(1,true);
     for(int i=0; i<n; i++) {  
        Rcpp::RawVector y=decompress(covresult[i]);
        //printf("LOOP %i NEW ",i);   
         //Rcpp::NumericVector y(covresult[i]); 
         int rangelen=y.size();
         // start from 1-based position (maybe Rcpp is 1-based)
         for(int k=1; k<= rangelen; k++){
          tempsum+=y[k];
         }

         if(tempsum>0){
          ALLzeros(0)=false;
          break;
         }
     }
     
     //ALLzeros(0)=allzeros;
     //return boolean. "true" if 
     return (ALLzeros); 
}












// [[Rcpp::export]]
NumericMatrix makeMatrixFrombaseCoverageCPP (List GRbaseCoverageOutput, int Nbins, int Snorm, IntegerVector keys) { 
      //find the size of the list of rawVectors (number of ranges)
      //the size of the list must be the same as the keys logicalVector
      int n = GRbaseCoverageOutput.size(); 
      //intermediate sum for each bin (very inside) 
      int res;
      Function decompress ("lz4_decompress_raw");
      //define final numeric matrix to give as result (nrows=n (number of ranges), ncols= Nbins(number of bins))
      Rcpp::NumericMatrix mat( n , Nbins );
      //vector to store length for each range, for final Snorm
      std::vector<int> lengths(n);

      //loop through each range in the list (i= index of the current range)
      for(int i=0; i<n; i++) { 

        //for this current range, extract the current raw vector (coverage vals...-> necessary?)
        //Rcpp::RawVector compressedrange(GRbaseCoverageOutput[i]); 

        //here decompress LZ4 each loop iteration: are we saving RAM?
        //I think this step is very inefficient, but still fast and acceptable
        Rcpp::RawVector currentrange=decompress(GRbaseCoverageOutput[i]);

        int m=currentrange.size();
        int currentkeychoice=keys[i];


        //HERE DIFFERENT KIND OF COMPRESSIONS for each range
        //1-byte compression 
        if (currentkeychoice==2){
          //divide m (size) by the number of bins: how many elements to be summed for each bin?
          int chunkblock=m/Nbins;
          // for each bin, sum elements in each bin. The remaining is considered 
          for(int k=0; k<Nbins; k++){
            res=0;
            for(int j=k*chunkblock; j<(k+1)*chunkblock; j++){     
              res+=currentrange[j]; 
            }   
            mat(i,k)=res;
          } 
          //here length range=m (the number of bp of the current range is exactly the number of bytes)
          //populate the array of lengths (parameter useful if Snorm=TRUE)
          lengths[i]=m;         
        }


        //nibble compression even/odd
        if (currentkeychoice==3 | currentkeychoice==4){
          //divide m (size) by the number of bins: how many elements to be summed for each bin?
          // if nibble compression, real chunk will be double of the number of bytes
          int chunkblock;
          int reminder;
          if(currentkeychoice==3){
            chunkblock=m*2/Nbins;
            //here length range=(m*2): number of bp is double of number of bytes
            //populate the array of lengths (parameter useful if Snorm=TRUE)
            lengths[i]=m*2; 
          }else{
            chunkblock=(m*2-1)/Nbins;
            reminder= (m*2-1)%Nbins;
            //here length range=(m*2)-1 : number of bp is double of number of bytes -1 (because extra byte)
            //populate the array of lengths (parameter useful if Snorm=TRUE)
            lengths[i]=(m*2)-1; 
          }
          
          int currenthalf=0;
          int counterblock=0;
          unsigned char num_byte;
          // for each bin, sum elements in each bin. The remaining is considered 
          for(int k=0; k<Nbins; k++){
            unsigned char missingpart;
            res=0;
            for(int j=k*chunkblock; j<(k+1)*chunkblock; j++){     
              if(currenthalf==0  ){
                num_byte= (currentrange[counterblock] & 0xF0) >>4 ;
              //else, take only last 4 bits (second nibble)
              }else{
                num_byte= currentrange[counterblock] & 0x0F ;     
              }
              res+=num_byte;

              currenthalf++;
              if (currenthalf==2){
                //reset, this byte expired
                counterblock++;
                currenthalf=0;
              }
            }

            //if 4 (odd) and nbin-1, add last byte of the currentrange (otherwise only most significant nibble is added (=0 by definition))
            //AND only if Nbins can be divided into chunks?? If remaining, do not do it
            if(currentkeychoice==4 & k==Nbins-1 &reminder==0){
              missingpart=currentrange[m-1] & 0x0F;
              res+=missingpart;
            }

            mat(i,k)=res;
          }
        }




        //2-byte compression 
        if (currentkeychoice==1){
          //cannot divide directly for Nbins, because the remaining part is distributed and we can have odd chunks
          int halfrange=m/2;
          int chunkblock=halfrange/Nbins;
          int doublechunkblock=chunkblock*2;
          for(int k=0; k<Nbins; k++){
            //here we are inside the bin
            int temp=0;
            //for each chunk (little array), sum everyhing 2 by 2 (2-byte compression, big endian)
            for(int j=k*doublechunkblock; j<(k+1)*doublechunkblock; j+=2){     
              temp+=currentrange[j]*256+currentrange[j+1];
            }   
            mat(i,k)=temp;
          }  
          //here length range=m/2 (n bp= half number of bytes)
          //populate the array of lengths (parameter useful if Snorm=TRUE)
          lengths[i]=m/2;          
        }
      


      }



      //if Snorm=TRUE (Snorm>0), divide each matrix value for the length of that range
      if(Snorm>0){
        for(int i=0; i<n; i++){
          for(int k=0; k<Nbins; k++){
            mat(i,k) = mat(i,k)/lengths[i];
          }
        }
     }
       
   return mat;  
}
