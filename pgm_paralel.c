#include "pgm_IO.h"
#include "pgm_IO.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi/mpi.h>

#define NRITER 8 // am folosit macro-ul pentru a defini lungimea vectorului cu valorile iteratiilor

int main(int argc, char ** argv)
{
    // Atribuim valorile din fisier lui M si N
    int M, N, nproc, MP, NP, rank;
    char *pgm_name_in = "image_640x480.pgm"; // dam numele fisierului ce contine datele de input
    
    pgm_size(pgm_name_in, &M, &N); // apelam functia de size pentru a prelua valorile lui M si N
    
    float masterdata[M][N]; // definim matricea masterdata care va retine fisierul grafic initial
    
    MPI_Init(&argc, &argv); // pornim mediul MPI
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // determinam rangul procesului si-l retinem in variabila rank
    MPI_Comm_size(MPI_COMM_WORLD, &nproc); // determinam numarul de procese asociate comunicatorului precizat si-l retinem in variabila nproc
    
    // PAS 1 -> Calculam MP si NP avand in vedere faptul ca imaginea este divizata de-a lungul liniilor de pixeli
    MP = M/nproc;
    NP = N;
    
    // PAS 2 -> In procesul master citim continutul fisierul grafic initial in masterdata folosind functia pgm_read
    if (rank == 0) 
    {
        pgm_read(pgm_name_in, masterdata, M, N);
    
    }
    
    /* Deoarece s-a mentionat ca toate procesele efectueaza secventele de cod ce urmeaza, n-am mai folosit o ramura else pentru ca ar fi trebuit sa 
     * duplicam cod.
    */
    
    // PAS 3 -> Declaram matricele necesare reconstructiei
    float pold[MP+2][NP+2], pnew[MP+2][NP+2], plim[MP+2][NP+2];
    float data[MP][NP];
    
    // PAS 4 
    /* Folosim MPI_Scatter pentru a distribui segmentele de date de dimensiune MP*NP                                                                                     
     * din masterdata catre toate procesele, adica le vom primi in matricea data. In cazul de fata, radacina este chiar procesul master, adica cel cu 
     * rang 0.
     */
    MPI_Scatter(masterdata, MP*NP, MPI_FLOAT, data, MP*NP, MPI_FLOAT, 0, MPI_COMM_WORLD); 
    
    // PAS 5 -> Initializam matricile  
    
    /* Am ales sa fac aceste for-uri inainte, deoarece este mai eficient pentru ca o sa initializam toate elementele matricelor, inclusiv zona de 
     * halo, si nu vom mai ajunge sa tratam cazuri mai speciale, precum cele pentru zona de halo, adica pentru valorile de pe liniile 0 si MP+1, 
     * respectiv coloanele 0 si NP+1.
     */
    
    int i,j;

    for (i=0;i<MP+2;i++)
        for (j=0;j<NP+2;j++)
            {
                pold[i][j] = pnew[i][j] = plim[i][j] = 255;
            }
    
    // Punem in plim valorile din data, adica valorile pixelilor imaginii initiale
    for (i=1;i<=MP;i++)
        for (j=1;j<=NP;j++)
            plim[i][j] = data[i-1][j-1];
    
    // PAS 5 -> Realizam algoritmul de reconstructie folosind formula data  
    int niter, idxiter = 0; // idxiter reprezinta un indice pe care-l folosim pentru a parcurge vectorul cu valorile iteratiilor
    int vect_iter[] = {20,200,400,600,800,1000,1200,1400}; /* vectorul cu valorile pe care trebuie sa le ia niter, adica numarul de iteratii pentru o                               
                                                              imagine*/
    
    while(idxiter<NRITER) /*parcurgem indicii pentru a prelua valori din vector pana cand depasim numarul de elemente din vector - 1, 
                            pentru ca am considerat ca incepem de la 0 cu indexarea*/
    {
        niter=vect_iter[idxiter]; // preluam numarul de iteratii din vector
        for (int iter=0;iter<niter;iter++)
        {
            /* Cele doua variabile, trimite_rank si recept_rank o sa retina rangul procesului catre care se trimite si rangul procesului de la care 
             * se receptioneaza. 
             * Cele doua if-uri sunt pentru a initializa cele doua ranguri cu un 'dummy' MPI_PROC_NULL in cazul in care rank-1 sau 
             * rank+1 reprezinta un rang inexistent, avand in vedere ca rank are valori intre 0 si nproc-1.
             */
            
            int trimite_rank = rank+1; 
            int recept_rank = rank-1; 
            
            if (trimite_rank>=nproc)
                trimite_rank = MPI_PROC_NULL;
            
            if (recept_rank<0)
                recept_rank = MPI_PROC_NULL;
            
            /* Utilizam MPI_Sendrecv pentru a trimite/a receptiona NP pixeli (lucram cu pixeli de margine si din regiunea de halo) catre/de la un proces.
             * De asemenea, pentru a nu include valorile pixelilor din zona de halo, am modificat bufferul data ca parametru, precizand si de la ce 
             * coloana sa inceapa. 
             */
            MPI_Sendrecv(&pold[MP][1], NP, MPI_FLOAT, trimite_rank, 0, &pold[0][1], NP, MPI_FLOAT, recept_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Am aplicat acelasi principiu, dar cu cateva valori inversate.
            trimite_rank=rank-1;
            recept_rank=rank+1;
            
            if (trimite_rank<0)
                trimite_rank = MPI_PROC_NULL;
            
            if (recept_rank>=nproc)
                recept_rank = MPI_PROC_NULL;
            
            MPI_Sendrecv(&pold[1][1], NP, MPI_FLOAT, trimite_rank, 0, &pold[MP+1][1], NP, MPI_FLOAT, recept_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (i=1;i<=MP;i++)
                for (j=1;j<=NP;j++)
                {   
                    pnew[i][j] = 0.25*(pold[i-1][j]+pold[i+1][j]+pold[i][j-1]+pold[i][j+1]-plim[i][j]); /*calculam valorile lui pnew pe baza 
                                                                                                          algoritmului de                                     
                                                                                                          reconstructie (un algortim de edge-detection inversat), adica reconstruim imaginea folosind structura de margini de zona*/
                }
            
            // copiem matricea pnew in pold fara valorile de halo 
            for (i=1;i<=MP;i++)
                for (j=1;j<=NP;j++)
                    pold[i][j]=pnew[i][j];
        }    
            
        // PAS 6 -> Copiem matricea pold in data fara valorile de halo
        for (i=1;i<=MP;i++)
            for (j=1;j<=NP;j++)
                data[i-1][j-1] = pold[i][j];
        
        // PAS 7 -> Adunam toate matricile data obtinute in matricea masterdata la procesul master cu MPI_Gather(parametrul root va fi 0)
        MPI_Gather(data, MP*NP, MPI_FLOAT, masterdata, MP*NP, MPI_FLOAT, 0, MPI_COMM_WORLD);    
            
        // PAS 8 -> Scriem fisierele pgm in cadrul procesului master
        if (rank == 0)
        {
            char pgm_name_out[128]; // stringul care va retine numele fisierelor pgm
            snprintf(pgm_name_out, 128, "image_640x480_paralel_%d.pgm", niter); /*adaug stringul cu numele fisierului de output in pgm_name_out,    
                                                                                  pentru a crea fisiere cu   
                                                                                  nume diferite, care sa contina si numarul iteratiilor*/
            pgm_write(pgm_name_out,masterdata,M,N); // apelam functia de write pentru a scrie matricea cu valorile pixelilor in fisierul respectiv
        }
        
        /* Deoarece nu dam numarul de iteratii ca parametru si astfel nu apelam programul de mai multe ori pentru a obtine imaginile ci o singura 
         * data va trebui sa reinitializam valorile lui pold la valoarea implicita 255. Pentru pnew nu este necesar, deoarece o sa fie suprascrise 
         * aceleasi valori, fiind matricea care este calculata, iar la pold este necesar pentru ca acesta este folosit in calculul lui pnew si o sa-
         * si pastreze valorile de la iteratia anterioara pentru calcul, ci nu 255.
         */
        for (i=0;i<MP+2;i++)
            for (j=0;j<NP+2;j++)
                {
                    pold[i][j] = 255;
                }
        
        idxiter++; // incrementam indicele pentru a trece la urmatorul numar de iteratii
    }
    
    // PAS 9 
    MPI_Finalize(); // Inchidem mediul MPI
    
    return 0; 
}
