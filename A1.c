//Includes
#include <stdio.h>
#include <gmp.h>
#include <math.h>
#include "mpi.h"

//Macro def for min func
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

int main(int argc, char **argv) {
   //MPI data init
   int numProc = 0;
   int pRank = 0;
   MPI_Status status;
   MPI_Request request;

   //Timer init
   double startTime = 0.0;
   double doneTime = 0.0;

   //Max num, should be either 1 billion or 1 trillion
   long long int totalLimit = 1000000000LL;

   //Start MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &pRank);
   MPI_Comm_size(MPI_COMM_WORLD, &numProc);

   //Long long int inits
   long long int previousPrime, firstPrime, secondPrime, currentResult, oldResult;
   oldResult = 0LL;

   //MPZ variable inits
   mpz_t currentPrime, lastPrime, nonePrime;
   mpz_inits(currentPrime, lastPrime, nonePrime);

   //Wait for each processor to finish inits
   MPI_Barrier(MPI_COMM_WORLD);

   //Start timer
   startTime = MPI_Wtime();

   //Figure out how much work each proc needs to do
   //Math, logic is pretty much straight from lecture
   long long int length;
   int startInd, remain;
   length = floor(totalLimit / numProc);
   if (pRank < totalLimit % numProc) {
      length++;
   }
   remain = MIN(pRank, totalLimit % numProc);
   startInd = pRank * length + remain;

   //Get prime from next proc to get overlapping primes completed
   mpz_set_ui(lastPrime, startInd + length);
   mpz_set_ui(nonePrime, startInd);
   mpz_nextprime(lastPrime, lastPrime);

   //First prime in length
   mpz_nextprime(currentPrime, nonePrime);
   previousPrime = mpz_get_ui(currentPrime);

   //Compare each successive prime set of differences against current greatest dif
   //Store values if new greatest dif
   while(mpz_cmp(lastPrime, currentPrime) >= 0 && mpz_get_ui(currentPrime) <= totalLimit) {
      if(mpz_cmp(lastPrime, currentPrime) >= 0 && mpz_get_ui(currentPrime) <= totalLimit) {
         currentResult = mpz_get_ui(currentPrime) - previousPrime;
         if(currentResult > oldResult) {
            oldResult = currentResult;
            firstPrime = previousPrime;
            secondPrime = mpz_get_ui(currentPrime);
         }
      previousPrime = mpz_get_ui(currentPrime);
      mpz_nextprime(currentPrime, currentPrime);
    }
   }

   if(pRank != 0) {
      //Send largest found gap and primes for those gaps to proc 0
      //Make array to send all the data
      long long int send[3] = {0LL, 0LL, 0LL};
      send[0] = firstPrime;
      send[1] = secondPrime;
      send[2] = oldResult;
      MPI_Send(send, 3, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
   }

   //First proc finding true largest gap from all procs
   if (pRank == 0) {
      //Set holding vars for best values
      long long int bestFirst, bestSecond, bestResult;
      bestFirst = 0LL;
      bestSecond = 0LL;
      bestResult = 0LL;
      long long int holder[3] = {0LL, 0LL, 0LL};

      //Iterate through all procs
      for (int src = 1; src < numProc; src++) { 

         //Get each proc's best candidate values
         MPI_Recv(holder, 3, MPI_LONG_LONG_INT, src, 0, MPI_COMM_WORLD, &status);

         //Decide if candidiate is best so far
         if (holder[2] > bestResult) {
            bestFirst = holder[0];
            bestSecond = holder[1];
            bestResult = holder[2];
         }
      }

      //Check for proc 0
      if (oldResult > bestResult) {
         bestFirst = firstPrime;
         bestSecond = secondPrime;
         bestResult = oldResult;
      }

      //Get end time
      doneTime = MPI_Wtime();

      //Print all results
      printf("\nTime taken: %.2lf seconds\n", doneTime - startTime);
      printf("Largest gap = %lld\n", bestResult);
      printf("First prime = %lld, second prime = %lld\n", bestFirst, bestSecond);
   }

   //End MPI
   MPI_Finalize();
   //Close/finish program
   return 0;
}