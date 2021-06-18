#include "pgm_IO.h"
#include "pgm_IO.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NRITER 8 // am folosit macro-ul pentru a defini lungimea vectorului cu valorile iteratiilor

int main(int argc, char ** argv)
{
    // PAS 1 -> atribuim valorile din fisier lui M si N
    int M, N;
    char *pgm_name_in = "image_640x480.pgm"; // dam numele fisierului ce contine datele de input
    
    pgm_size(pgm_name_in, &M, &N); // apelam functia de size pentru a prelua valorile lui M si N
    
    // PAS 2 -> citim din fisierul declarat anterior matricea
    
    // Declaram matricele necesare reconstructiei
    float pold[M+2][N+2], pnew[M+2][N+2], plim[M+2][N+2];
    float data[M][N];
    
    pgm_read(pgm_name_in, data, M, N); // apelam functia de read pentru a citi continutul fisierului in matricea data
    
    /* PAS 4 -> Am ales sa fac acest pas inainte de 3, deoarece este mai eficient pentru ca o sa initializam toate elementele matricelor, inclusiv zona de halo, 
     * si nu vom mai ajunge sa tratam cazuri mai speciale, precum cele pentru zona de halo, adica pentru valorile de pe liniile 0 si MP+1, respectiv coloanele 0 si NP+1.*/
    int i,j;

    for (i=0;i<M+2;i++)
        for (j=0;j<N+2;j++)
            {
                pold[i][j] = pnew[i][j] = plim[i][j] = 255;
            }
    
    // PAS 3 -> Punem in plim valorile din data, adica valorile pixelilor imaginii initiale
    for (i=1;i<=M;i++)
        for (j=1;j<=N;j++)
            plim[i][j] = data[i-1][j-1];
    
    // PAS 5 -> Realizam algoritmul de reconstructie folosind formula data
    int niter, idxiter = 0; // idxiter reprezinta un indice pe care-l folosim pentru a parcurge vectorul cu valorile iteratiilor
    int vect_iter[] = {20,200,400,600,800,1000,1200,1400}; // vectorul cu valorile pe care trebuie sa le ia niter, adica numarul de iteratii pentru o imagine
    
    while(idxiter<NRITER) /*parcurgem indicii pentru a prelua valori din vector pana cand depasim numarul de elemente din vector-1, 
                            pentru ca am considerat ca incepem de la 0 cu indexarea*/
    {
        niter=vect_iter[idxiter]; // preluam numarul de iteratii din vector
        for (int iter=0;iter<niter;iter++)
        {
            for (i=1;i<=M;i++)
                for (j=1;j<=N;j++)
                    pnew[i][j] = 0.25*(pold[i-1][j]+pold[i+1][j]+pold[i][j-1]+pold[i][j+1]-plim[i][j]); /*calculam valorile lui pnew pe baza algoritmului de reconstructie      
                                                                                                          (un algortim de edge-detection inversat),
                                                                                                          adica reconstruim imaginea folosind structura de margini de zona */
        
            // copiem matricea pnew in pold fara valorile de halo    
            for (i=1;i<=M;i++)
                for (j=1;j<=N;j++)
                    pold[i][j] = pnew[i][j];
        }    
            
        // PAS 7 -> copiem matricea pold in data fara valorile de halo
        for (i=1;i<=M;i++)
            for (j=1;j<=N;j++)
                data[i-1][j-1] = pold[i][j];
        
        // PAS 8 -> scriem fisierele pgm
        char pgm_name_out[128]; // stringul care va retine numele fisierelor pgm
        snprintf(pgm_name_out, 128, "image_640x480_serial_%d.pgm", niter); /*adaug stringul cu numele fisierului de output in pgm_name_out, pentru a crea fisiere cu nume 
                                                                             diferite, care sa contina si numarul iteratiilor*/
        pgm_write(pgm_name_out,data,M,N); // apelam functia de write pentru a scrie matricea cu valorile pixelilor in fisierul respectiv
        
        /* Deoarece nu dam numarul de iteratii ca parametru si astfel nu apelam programul de mai multe ori pentru a obtine imaginile, ci o singura data
         * va trebui sa reinitializam valorile lui pold la valoarea implicita 255. Pentru pnew nu este necesar, deoarece o sa fie suprascrise aceleasi valori, fiind 
         * matricea care este calculata, iar la pold este necesar pentru ca acesta este folosit in calculul lui pnew si o sa-si pastreze valorile de la iteratia 
         * anterioara pentru calcul, ci nu 255.
         */
        for (i=0;i<M+2;i++)
            for (j=0;j<N+2;j++)
                {
                    pold[i][j] = 255; 
                }
        
        idxiter++; // incrementam indicele pentru a trece la urmatorul numar de iteratii
    }
    
    return 0;
}
