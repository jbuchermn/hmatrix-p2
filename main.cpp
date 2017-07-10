#include "hMatrixP2.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <ctime>

#include <Eigen/Core>
#include <Eigen/Dense>

void get_matrix(std::complex<double>* mat, int dim, HMatrixP2* hmat){
    std::vector<std::complex<double> > right(dim);
    std::vector<std::complex<double> > left(dim);

    for(int j=0,ij=0;j<dim;j++){
        for(unsigned int k=0; k<dim; k++) left[k]=0;
        for(unsigned int k=0; k<dim; k++) right[k]=0;
        right[j]=1.;
        hmat->apply(left.data(), right.data());

        for(int k=0;k<dim;k++,ij++)
            mat[ij]=left[k];
    }
}

void mult(const std::complex<double>*  matrix, std::complex<double>*  left_vec, const std::complex<double>*  right_vec, unsigned int dim){
    for(unsigned int i=0; i<dim; i++){
        for(unsigned int j=0; j<dim; j++){
            left_vec[j]+=matrix[j+dim*i]*right_vec[i];
        }
    }
}


int main(int argc, char* argv[]){
    int dim = 256;
    int runs= 1000;
    if(argc>1){
        dim=std::stoi(argv[1]);
    }

    std::cout<<"H-Matrix test, dim="<<dim<<std::endl;

    std::vector<std::complex<double> > matrix(dim*dim);
    std::vector<std::complex<double> > left_vec(dim);
    std::vector<std::complex<double> > right_vec(dim);

    for(unsigned int i=0; i<dim*dim; i++) matrix[i]    = std::complex<double>(1./(i+1), i);
    for(unsigned int i=0; i<dim;     i++) right_vec[i] = std::complex<double>(i, 1./(i+1));
    
    clock_t t;
    double ns;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(unsigned int i=0; i<dim; i++) left_vec[i]=0.;
    t = clock();
    for(unsigned int i=0; i<runs; i++) mult(matrix.data(),left_vec.data(),right_vec.data(),dim);
    t = clock() - t;
    ns = 1.e9*double(t)/CLOCKS_PER_SEC;
    std::cout<<"Regular application:  "<<ns/runs<<"ns, "<<(ns/runs/dim/dim)<<"ns, "<<dim*dim<<std::endl;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(unsigned int i=0; i<dim; i++) left_vec[i]=0.;
    t = clock();
    for(unsigned int i=0; i<runs; i++) Eigen::Map<Eigen::VectorXcd>(left_vec.data(), dim)+=Eigen::Map<Eigen::MatrixXcd>(matrix.data(), dim, dim)*Eigen::Map<Eigen::VectorXcd>(right_vec.data(),dim);
    t = clock() - t;
    ns = 1.e9*double(t)/CLOCKS_PER_SEC;
    std::cout<<"Eigen application:    "<<ns/runs<<"ns, "<<(ns/runs/dim/dim)<<"ns, "<<dim*dim<<std::endl;
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    HMatrixP2 hmat(matrix.data(), HMatrixP2::IndexPartitioning(dim), HMatrixP2::IndexPartitioning(dim));

    std::vector<std::complex<double> > check(dim*dim);
    get_matrix(check.data(), dim, &hmat);
   
    double frob_error=0.;
    double frob_mat=0.;
    for(unsigned int i=0; i<dim; i++){
        for(unsigned int j=0; j<dim; j++){
            frob_mat   += std::pow(std::abs(matrix[j+dim*i]               ),2);
            frob_error += std::pow(std::abs(check[j+dim*i]-matrix[j+dim*i]),2);
        }
    }
    std::cout<<std:endl<<"HMatrix truncation error rel/abs: "<<std::sqrt(frob_error/frob_mat)<<"/"<<std::sqrt(frob_error)<<std::endl;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(unsigned int i=0; i<dim; i++) left_vec[i]=0.;
    t = clock();
    for(unsigned int i=0; i<runs; i++) hmat.apply(left_vec.data(),right_vec.data());
    t = clock() - t;
    ns = 1.e9*double(t)/CLOCKS_PER_SEC;
    std::cout<<"HMatrix application:  "<<ns/runs<<"ns, "<<(ns/runs/hmat.applyCount())<<"ns, "<<hmat.applyCount()<<std::endl;

}
