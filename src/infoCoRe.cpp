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
  //arma::SpMat<double> L(N,N); L=D-Adj;
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
void test2( const arma::SpMat<double>& Adj, Rcpp::IntegerVector weighted=0){

  // Refs: 1) https://stackoverflow.com/questions/67189074/hermitian-adjacency-matrix-of-digraph
  // Refs: 2) https://stackoverflow.com/questions/67189074/hermitian-adjacency-matrix-of-digraph
  // Refs: 3) https://www.stats.ox.ac.uk/~cucuring/Hermitian_Clustering_AISTATS.pdf  

  arma::uword i,j,ii,jj,N,option_we;

  option_we = weighted[0];
  
  Adj.brief_print("Adj:");
  
  N = Adj.n_rows;  
  
  arma::SpMat<std::complex<double>> H; H.zeros(N,N);

  if( option_we == 0 ){
    //unweighted    
  
    for(ii=0; ii<N; ii++){
      const arma::SpSubview_row<double> rindx = Adj.row(ii);
      const arma::uvec cindx = find(rindx);
      for(jj=0; jj<cindx.n_elem; jj++){

        // Define the Adj. matrix rows and column indices.
        i = ii;
        j = cindx(jj);

        // Find edge weights
        double we_ij=Adj.at(i,j);
        double we_ji=Adj.at(j,i);

        // Build Hermitain matrix
        if( we_ij > 0 && we_ji > 0 ){
          // edge in both directions
          std::complex<double> bidir(1,0);
          H.at(i,j) = bidir;
        } else {
          if( we_ij > 0 && we_ji == 0 ){
            // out-going edge
            std::complex<double> cWe_ij(0,1);
            H.at(i,j) = cWe_ij; H.at(j,i) = std::conj(cWe_ij);
          } else { 
            // in-coming edge
            std::complex<double> cWe_ji(0,1);        
            H.at(j,i) = cWe_ji; H.at(i,j) = std::conj(cWe_ji);
          }
        }
      }
    }

  } else {
    //weighted
    
     for(ii=0; ii<N; ii++){
      const arma::SpSubview_row<double> rindx = Adj.row(ii);
      const arma::uvec cindx = find(rindx);
      for(jj=0; jj<cindx.n_elem; jj++){

        // Define the Adj. matrix rows and column indices.
        i = ii;
        j = cindx(jj);

        // Find edge weights
        double we_ij=Adj.at(i,j);
        double we_ji=Adj.at(j,i);

        // Build Hermitain matrix
        if( we_ij > 0 && we_ji > 0 ){
          // edge in both directions
          std::complex<double> bidir((we_ij-we_ji),0);
          H.at(i,j) = bidir;
        } else {
          if( we_ij > 0 && we_ji == 0 ){
            // out-going edge
            std::complex<double> cWe_ij(0,we_ij);
            H.at(i,j) = cWe_ij; H.at(j,i) = std::conj(cWe_ij);
          } else { 
            // in-coming edge
            std::complex<double> cWe_ji(0,we_ji);        
            H.at(j,i) = cWe_ji; H.at(i,j) = std::conj(cWe_ji);
          }
        }
      }
    }
     
  }


  //Check H is Hermitain
  cout << "Is H hermitian:  " << H.is_hermitian() << endl;
  cout << "Is H symmetric:  " << H.is_symmetric() << endl;
  //cout << "Is H her & symm: " << H.is_sympd() << endl;
   
  H.brief_print("H:");

  arma::SpMat<double> d_in  = sum(Adj,0);
  arma::SpMat<double> d_out = sum(Adj,1);
  arma::SpMat<double> d_all = d_in.t() + d_out;
  arma::SpMat<double> D_r(N,N); D_r.diag() = d_all;
  arma::SpMat<double> D_i(N,N); D_i.zeros();
  
  arma::SpMat<std::complex<double>> L(D_r,D_i);

  L.brief_print("L (before):");
  
  // Laplacian
  L = L-H;
  //arma::SpMat<std::complex<double>> L(N,N); L=L.diag(D.diag())-H;

  L.brief_print("L (after):");
  
  cout << "Is L hermitian:  " << L.is_hermitian() << endl;
  cout << "Is L symmetric:  " << L.is_symmetric() << endl;
  //cout << "Is L her & symm: " << L.is_sympd() << endl;
  
  L.brief_print("L:");

  d_in.transform(  [](double val) { return ( (val==0 ? val : pow(val,-0.5)) ); });
  d_out.transform( [](double val) { return ( (val==0 ? val : pow(val,-0.5)) ); });

  // Normalised Laplacian  
  arma::SpMat<double> D_in_r(N,N);  D_in_r.diag()  = d_in;
  arma::SpMat<double> D_out_r(N,N); D_out_r.diag() = d_out;

  arma::SpMat<std::complex<double>> I; I.eye(N,N);
  arma::SpMat<std::complex<double>> D_in (D_in_r, D_i);
  arma::SpMat<std::complex<double>> D_out(D_out_r,D_i);
  
  arma::SpMat<std::complex<double>> L_hat = D_in * (H+I) * D_out;

  L_hat.brief_print("L_hat:");

  cout << "Is L_hat hermitian:  " << L_hat.is_hermitian() << endl;
  cout << "Is L_hat symmetric:  " << L_hat.is_symmetric() << endl;
  //cout << "Is L_hat her & symm: " << L_hat.is_sympd() << endl;

  L_hat.brief_print("L_hat:");
}
