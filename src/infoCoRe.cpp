#include "Headers.h"

//----------------
//NOTE:
// When editing this file, i.e. adding or removing functions, need to carry out the following tasks:
// 1) edit file NAMESPACE, adding/removing the function name in export(...)
// Next run following in R, inside the package:
// 2) $> cd rblock/
// 3) $> R
// 4)  > library(Rcpp)
// 5)  > Rcpp::compileAttributes()
// running 5) will create new RcppExports.R in R/ and
// new RcppExports.cpp in src/
// If need to run (to build Makevars file in src):
// 1) $> cd rblock/
// 2) $> autoconfig
// 3) $> ./configure
// BUILD PACKAGE:
// 1) $> R CMD build rblock/
// INSTALL PACKAGE (note 'INSTALL' should be in capital letters):
// 1) $> R CMD INSTALL rblock_1.0.tar.gz
//----------------


//Ref Papers:
// 1) P. Villegas et al. Laplacian paths in complex networks: Information core emerges from entropic transitions, PHYSICAL REVIEW RESEARCH 4, 033196 (2022).
// doi: DOI: 10.1103/PhysRevResearch.4.033196
// website: https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.4.033196

// [[Rcpp::export]]
void test( const arma::SpMat<double>& Adj){
//void test( const arma::Mat<double>& Adj){

  // TEST:
  // > library(infoCoRe)
  // > M1<-matrix(rpois(36,5),nrow=6)
  // > infoCoRe::test(Adj=M1)  
  //Adj.print();
  
  uint_fast32_t i,N;

  N = Adj.n_rows;

  /*
  arma::rowvec d_in  = sum(Adj,0);
  arma::colvec d_out = sum(Adj,1);
  arma::colvec d_all = d_in + d_out;
  */
  arma::SpMat<double> d_in  = sum(Adj,0);
  arma::SpMat<double> d_out = sum(Adj,1);
  arma::SpMat<double> d_all = d_in.t() + d_out;

  arma::SpMat<double> D(N,N); D.diag() = d_all;

  d_in.brief_print("d.in:");
  d_out.brief_print("d.out:");

  D.brief_print("D:");

  // Laplacian
  arma::SpMat<double> L(N,N); L=D-Adj;

  L.brief_print("L:");
  //cout << d_all.n_rows << ", " << d_all.n_cols << endl;
  //d_all.print();

  d_in.transform(  [](double val) { return ( (val==0 ? val : pow(val,-0.5)) ); });
  d_out.transform( [](double val) { return ( (val==0 ? val : pow(val,-0.5)) ); });

  d_in.brief_print("d.in:");
  d_out.brief_print("d.out:");

  // Normalised Laplacian
  arma::SpMat<double> I; I.eye(N,N);
  arma::SpMat<double> D_in(N,N);  D_in.diag()  = d_in;
  arma::SpMat<double> D_out(N,N); D_out.diag() = d_out;

  I.brief_print("I:");
  D_in.brief_print("D.in:");
  D_out.brief_print("D.out:");
  
  arma::SpMat<double> L_hat = D_in * (Adj+I) * D_out;

  L_hat.brief_print("L.hat:");
  
  cout << "\n";

  /*
  cout << "TEST find indices of non-zero ele in Adj: " << endl;

  arma::uword ii, jj, NR;
  NR = Adj.n_cols;

  vector<EDGE> edges;
  
  for(ii=0; ii<NR; ii++){
    const arma::SpSubview_col<double> cc = Adj.col(ii);
    const arma::uvec rindx = find(cc);
    //p2.brief_print();
    for(jj=0; jj<rindx.n_elem; jj++){
      edges.push_back( EDGE(ii,rindx(jj), Adj(ii,rindx(jj))) );
    }
  }

  cout << "edges size: " << edges.size() << endl;

  for(i=0; i<20; i++){
    cout << "[" << std::get<0>(edges[i]) << "," << std::get<1>(edges[i]) << "]" << endl;
  }
  
  //free edges memory vector
  std::vector<EDGE>().swap(edges);
  */
 
  
}

// [[Rcpp::export]]
void test2( const arma::SpMat<double>& Adj){

  // Refs: 1) https://www.stats.ox.ac.uk/~cucuring/Hermitian_Clustering_AISTATS.pdf  
 
  Adj.brief_print("Adj:");

  arma::uword i,j,ii,jj,N;
  N = Adj.n_rows;

  arma::SpMat<std::complex<double>> H; H.zeros(N,N);
  
  for(ii=0; ii<N; ii++){
    const arma::SpSubview_row<double> rindx = Adj.row(ii);
    const arma::uvec cindx = find(rindx);
    for(jj=0; jj<cindx.n_elem; jj++){

      // Define the Adj. matrix rows and column indices.
      i = ii;
      j = cindx(jj);

      //std::complex<double> we(0,Adj.at(i,j));
      //H.at(i,j) = we;
      
      
      // Deine the edge weight as a complex number
      std::complex<double> we(0,Adj.at(i,j));
      
      if( i <= j ){
        // if out-going edges, or edges in both directions

        // edges in both directions
        if( i == j ){ H.at(i,j) = we*std::conj(we); }
        // out-going edges
        else        { H.at(i,j) = we; H.at(j,i) = std::conj(we); }
        // in-coming edges
      } else        { H.at(j,i) = we; H.at(i,j) = std::conj(we); }
      
    }
  }


  //Check H is Hermitain
  cout << "Is H hermitian: " << H.is_hermitian() << endl;

    H.brief_print("H:");

}
