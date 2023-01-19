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
    
  cout << d_all.n_rows << ", " << d_all.n_cols << endl;
  //d_all.print();

  cout << "\n";

  
  //arma::Mat<double> D(N,N);
  //for(i=0; i<N; i++){
  //  Adj.row(i);
  //}

  
  /*
  // pointer to matrix object on disk
  XPtr<FBM_RW> xSmat = Smat["address_rw"];
  BMAcc_RW<double> Sr(xSmat); 

  //if(model==nullptr){ cout << "create model first." << endl; return; }
  
  arma::mat Y  = FBM2arma( Pmat ); 

  
  //--- Setup OpenMP
  //model->setOpenMP( nCORES[0] );

  //select the dataset
  int dataset = Dataset[0];

  uint_fast32_t i,j,k;

  double val;
  
  //size of Smat
  uint_fast32_t N = Sr.nrow();
  
  //set matrix (with diagonal) size
  uint_fast32_t K = N*(N+1)/2;

  //get a mapping from 2d matrix (with diagonal) to 1d array
  if( model->map2Arr.empty() ){
    //vector<tripleInt> map2Arr(K);
    model->setMapSpace(K); 
    model->mapMatrixDiag2Array( N, model->map2Arr );
  }

  //calculate t(Sr) * Y * Sr => t(CSFcsr) * t(Y) = res, res2 = res * CSFcsc  
  vector<tupleCOO> Srmap, tSrmap;
 
  CSF* Srcsr  = new CSF(); //store Sr1 transpose in csr format
  CSF* Srcsc  = new CSF(); //store Sr1           in csc format
  
  model->mapSrmatrix( Sr, N, Srmap, model->map2Arr, true );//read Sr in column-major format

  tSrVec( Srmap, tSrmap, N );//transpose

  cout << "Y" << endl;
  cout << Y << endl;

  cout << "tSrmap" << endl;
  for(k=0; k<tSrmap.size(); k++){
    i = std::get<0>(tSrmap[k]);
    j = std::get<1>(tSrmap[k]);
    val = std::get<2>(tSrmap[k]);
    cout << "(" << i << "," << j << ") = " << val << endl;
  }

    
  model->mappingCOO2CSF( tSrmap, Srcsr, N, false );//store Sr in CSR format
  model->mappingCOO2CSF(  Srmap, Srcsc, N, true  );//store Sr in CSC format
 
  cout << "Srcsr" << endl;
  for(i=0; i<Srcsr->nz;i++){
    cout << "[" << i << "] val: " << Srcsr->val[i] << ", col: " << Srcsr->col[i] << endl; 
  }

  for(i=0; i<(Srcsr->nrows+1);i++){
    cout << "[" << i << "] row: " << Srcsr->row[i] << endl; 
  }
  cout << "---" << endl;


  arma::mat res,res2;
  model->csr_dense_tcrossprod( Srcsr,   Y, res   );
  model->dense_csc_prod      ( res, Srcsc, res2  );

  cout << "tSr * Y * Sr" << endl;
  cout << res2 << endl;


  arma::SpMat<double> X; X.zeros(N,N);
  SrVec2ArmaD( X, Srmap, N );

  arma::mat res3;
  cout << "tSr * Y * Sr arma " << endl;
  res3 = X.t() * Y * X;
  cout << res3 << endl;

  bool same = approx_equal(res2, res3, "reldiff", 0.1);

  if( same ){
    cout << "matices are the same" << endl;
  } else {
    cout << "matices are different" << endl;
  }
   

  //free space
  model->freeMapSpace();
  std::vector<tupleCOO>().swap(Srmap);
  std::vector<tupleCOO>().swap(tSrmap); 
  if(Srcsr){ delete Srcsr; }
  if(Srcsc){ delete Srcsc; }
  */
  
}

