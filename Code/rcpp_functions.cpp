#include <Rcpp.h>
#include <numeric>
#include <iostream>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
// MVP with Kronecker Matrix: Ax, with A = A_list[[1]] otimes ... otimes A_list[[P]] with P leq 4
NumericVector MVP_kronecker_rcpp(const List& A_list, const NumericVector& x){
  
  int P=A_list.size();
  std::vector<int> nrow_p(P);
  std::vector<int> ncol_p(P);
  int max_abs=1;
  int length_res=1;
  
  for(int p=P-1; p>=0; p--){
    const NumericMatrix& A_p=A_list[p];
    nrow_p[p]=A_p.nrow();
    ncol_p[p]=A_p.ncol();
    if(nrow_p[p]>ncol_p[p]){
      max_abs*=nrow_p[p];
    } else{
      max_abs*=ncol_p[p];
    }
    length_res*=nrow_p[p];
  }
  
  std::vector<double> v1(max_abs),v2(max_abs);
  
  for(int i=0; i<x.size(); i++){
    v1[i] = x[i];
  }
  
  for(int p=P-1; p>=0; p--){
    
    const NumericMatrix& A_p=A_list[p];
    int left=1;
    int right=1;
    for(int i=0; i<p; i++){
      left *= ncol_p[i];
    }
    for(int i=P-1; i>p; i--){
      right *= nrow_p[i];
    }
    int nrow=nrow_p[p];
    int ncol=ncol_p[p];
    int base_col=0;
    int base_row=0;
    
    v1.swap(v2);
    
    for(int s=0; s<left; s++){
      for(int j=0; j<right; j++){
        int index_col=base_col+j;
        NumericVector z_in(ncol);
        for(int l=0; l<ncol; l++){
          z_in[l]=v2[index_col];
          index_col += right;
        }
        NumericVector z_out(nrow);
        for(int r=0; r<nrow; r++){
          for(int c=0; c<ncol; c++){ 
            z_out[r] += A_p(r,c)*z_in[c];
          }
        }
        int index_row = base_row+j;
        for(int l=0; l<nrow; l++){
          v1[index_row] = z_out[l];
          index_row += right;
        }
      }
      base_col += (right*ncol);
      base_row += (right*nrow);
    }
    
  }
  
  NumericVector res(length_res);
  for(int i=0; i<res.size(); i++){
    res[i] = v1[i];
  }
  return res;
}


// [[Rcpp::export]]
// Diagonal of Kronecker Matrix: diag_A, A = A_list[[1]] otimes ... otimes A_list[[P]] with P leq 44
NumericVector diag_kronecker_rcpp(const List& A_list){
  
  int P = A_list.size();
  int K=1;
  std::vector<int> K_p(P);
  for(int p=0; p<P; p++){
    const NumericMatrix& A_p = A_list[p];
    K_p[p]=A_p.ncol();
    K*=K_p[p];
  }
  NumericVector diag(K);
  
  if(P==1){
    int j=0;
    const NumericMatrix& A_0 = A_list[0];
    for(int m_0=0; m_0<K_p[0]; m_0++){
      diag[j++] = A_0(m_0,m_0);
    }
    return diag;
  }
  else if(P==2){
    int j=0;
    const NumericMatrix& A_0 = A_list[0];
    const NumericMatrix& A_1 = A_list[1];
    for(int m_0=0; m_0<K_p[0]; m_0++){
      const double a0_m0_m0 = A_0(m_0,m_0);
      for(int m_1=0; m_1<K_p[1]; m_1++){
        diag[j++] = a0_m0_m0*A_1(m_1,m_1);
      }
    }
    return diag;
  }
  else if(P==3){
    int j=0;
    const NumericMatrix& A_0 = A_list[0];
    const NumericMatrix& A_1 = A_list[1];
    const NumericMatrix& A_2 = A_list[2];
    for(int m_0=0; m_0<K_p[0]; m_0++){
      const double a0_m0_m0 = A_0(m_0,m_0);
      for(int m_1=0; m_1<K_p[1]; m_1++){
        const double a1_m1_m1 = A_1(m_1,m_1);
        for(int m_2=0; m_2<K_p[2]; m_2++){
          diag[j++] = a0_m0_m0*a1_m1_m1*A_2(m_2,m_2);
        }
      }
    }
    return diag;
  }
  else if(P==4){
    int j=0;
    const NumericMatrix& A_0 = A_list[0];
    const NumericMatrix& A_1 = A_list[1];
    const NumericMatrix& A_2 = A_list[2];
    const NumericMatrix& A_3 = A_list[3];
    for(int m_0=0; m_0<K_p[0]; m_0++){
      const double a0_m0_m0 = A_0(m_0,m_0);
      for(int m_1=0; m_1<K_p[1]; m_1++){
        const double a1_m1_m1 = A_1(m_1,m_1);
        for(int m_2=0; m_2<K_p[2]; m_2++){
          const double a2_m2_m2 = A_2(m_2,m_2);
          for(int m_3=0; m_3<K_p[3]; m_3++){
            diag[j++] = a0_m0_m0*a1_m1_m1*a2_m2_m2*A_3(m_3,m_3);
          }
        }
      }
    }
    return diag;
  }
  else{
    std::cout<<"P too large"<<std::endl;
    return 0;
  }
  
}


 
