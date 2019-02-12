#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
Nagy Péter, M07ILF
2017.04.23.
*/

const double threshold=0.000000000000000001;

struct StrMtx{
    int width;
    int height;
    double* Mtx;
    int size;
};

struct IntStrMtx{
    int width;
    int height;
    int* Mtx;
    int size;
};

struct Element{
    int RowPlace;
    int ColPlace;
    double value;
};

struct StrElimMtx{
    struct StrMtx Left;
    struct StrMtx Right;
};

void Usage(){
    printf("Tobbvaltozos polinomillesztes: ./a.out (valtozok szama) (polinom fokszama) adatsor [kimeneti file-nevek (opcionalisan)] \n");
    exit(0);
}

struct StrMtx ReadMatrix(char* filename, long size){
    FILE *inMatrixStream;
    double *matrix;
    long msize=20*size*size;
    matrix=(double *)malloc(msize*sizeof(double));
    inMatrixStream=fopen(filename, "r");
    long counter=0;
    long InRowCounter=1;
    long rowcounter=1;
    while(1){
    if(feof(inMatrixStream)){
        break;
    }
        fscanf(inMatrixStream, "%lf", &matrix[counter]);
        if(counter==msize){
            msize=2*msize;
            matrix=(double *) realloc(matrix, msize*sizeof(double));
        }
        counter++;
    }
    struct StrMtx Result;
    Result.height=counter/size;
    Result.width=size;
    Result.Mtx=matrix;

    return Result;
}

void PrintMatrix(struct StrMtx inputmatrix, char* filename){
    FILE *outstream;
    outstream=fopen(filename, "w");
    int i,j;
    for(i=0;i<inputmatrix.height;i++){
        for(j=0;j<inputmatrix.width;j++){
            fprintf(outstream, "%f ", inputmatrix.Mtx[inputmatrix.width*i+j]);
        }
        fprintf(outstream, "\n");
    }
    fclose(outstream);
}

void PrintIntMatrix(struct IntStrMtx inputmatrix, char* filename){
    FILE *outstream;
    outstream=fopen(filename, "w");
    int i,j;
    for(i=0;i<inputmatrix.height;i++){
        for(j=0;j<inputmatrix.width;j++){
            fprintf(outstream, "%d", inputmatrix.Mtx[inputmatrix.width*i+j]);
        }
        fprintf(outstream, "\n");
    }
    fclose(outstream);
}

struct StrMtx TransposeMatrix(struct StrMtx matrix){
    int i, j;
    struct StrMtx transposed;
    int width=matrix.width;
    int height=matrix.height;
    double* transp;
    transp=(double *)malloc(width*height*sizeof(double));
    for(i=0; i<height; i++){
        for(j=0;j<width; j++){
            transp[j*height+i]=matrix.Mtx[i*width+j];
        }
    }
    transposed.Mtx=transp;
    transposed.width=height;
    transposed.height=width;
    return transposed;

}

struct StrMtx Szorzas(struct StrMtx matrixleft, struct StrMtx matrixright){
    int i, j, k;
    if(matrixleft.width!= matrixright.height){
        printf("A matrixok  nem kompatibilisek. \n");
        exit(0);
    }
    double *product=malloc(matrixleft.height*matrixright.width*sizeof(double));
    for(i=0;i<matrixleft.height;i++){
        for(j=0;j<matrixright.width;j++){
            double sum=0.0;
            for(k=0;k<matrixleft.width;k++){
                sum=sum+matrixleft.Mtx[i*matrixleft.width+k]*matrixright.Mtx[k*matrixright.width+j];
            }
            product[i*matrixright.width+j]=sum;
        }
    }
    struct StrMtx Result;
    Result.width=matrixright.width;
    Result.height=matrixleft.height;
    Result.Mtx=product;
    return Result;
}

struct IntStrMtx GenerateCoeffs(int NumberOfVariables, int Degrees){
    int MaxCount=pow(Degrees+1, NumberOfVariables)-1;
    struct IntStrMtx Result;
    Result.width=NumberOfVariables;
    int height=0;
    int *resultmatrix;
    resultmatrix=(int *)malloc(NumberOfVariables*MaxCount*sizeof(int));
    int i;


    for(i=0;i<MaxCount;i++){
        int *dummy;
        dummy=(int *)malloc(NumberOfVariables*sizeof(int));
        dummy[0]=i;
        int j;
        for(j=0;j<NumberOfVariables;j++){
            if(j+1<NumberOfVariables){
                dummy[j+1]=0;
                while(dummy[j]>(Degrees)){
                    dummy[j]=dummy[j]-(Degrees+1);
                    dummy[j+1]=dummy[j+1]+1;
                }
            }
            else{
                while(dummy[j]>(Degrees)){
                    dummy[j]=dummy[j]-(Degrees+1);
                }
            }
        }
        int sum=0;
        int k;
        for(k=0;k<NumberOfVariables;k++){
            sum=sum+dummy[k];
        }
        if(sum<=Degrees){
            int p;
            for(p=0;p<NumberOfVariables;p++){
                resultmatrix[height*NumberOfVariables+p]=dummy[p];
            }
            height=height+1;
        }
    }
    int *ResMatrix;
    ResMatrix=(int *)malloc(height*NumberOfVariables*sizeof(int));
    for(i=0;i<height*NumberOfVariables;i++){
        ResMatrix[i]=resultmatrix[i];
    }
    free(resultmatrix);
    Result.height=height;
    Result.Mtx=ResMatrix;
    return Result;
}

struct StrMtx GenerateVector(struct StrMtx Initial){
    int n=Initial.width-2;
    int N=Initial.height;
    int i;
    double* vector;
    vector=(double *)malloc(N*sizeof(double));
    for(i=0; i<N; i++){
        vector[i]=Initial.Mtx[i*(n+2)+n]/Initial.Mtx[i*(n+2)+n+1];
    }
    struct StrMtx Result;
    Result.width=1;
    Result.height=N;
    Result.Mtx=vector;
    return Result;
}

struct StrMtx GenerateErrorVector(struct StrMtx Initial){
    int n=Initial.width-2;
    int N=Initial.height;
    int i;
    double* vector;
    vector=(double *)malloc(N*sizeof(double));
    for(i=0; i<N; i++){
        vector[i]=Initial.Mtx[i*(n+2)+n+1];
    }
    struct StrMtx Result;
    Result.width=1;
    Result.height=N;
    Result.Mtx=vector;
    return Result;
}

struct StrMtx GenerateTervmtx(struct StrMtx Initial, struct IntStrMtx BaseFunctions){
    int nvar=Initial.width-2;
    int nbase=BaseFunctions.height;
    int N=Initial.height;
    int width=nbase;
    double* tervmatrix;
    tervmatrix=(double *)malloc(nbase*N*sizeof(double));
    int i,j, k,l;
    for(k=0; k<N; k++){
        double sigm=1./Initial.Mtx[k*(nvar+2)+nvar+1];
        for(i=0;i<nbase;i++){
            double estimate=1.0;
            for(l=0;l<nvar;l++){
                if(BaseFunctions.Mtx[i*nvar+l]!=0){
                    estimate=estimate*pow(Initial.Mtx[k*(nvar+2)+l], BaseFunctions.Mtx[i*nvar+l]);
                }
            }
            tervmatrix[k*(width)+i]=sigm*estimate;
        }
    }
    struct StrMtx Result;
    Result.height=N;
    Result.width=width;
    Result.Mtx=tervmatrix;
    return Result;
}

struct StrMtx IdN(int n){
    struct StrMtx idn;
    double *idnmatrix;
    int i;
    int j;
    idnmatrix= (double*) malloc(sizeof(double)*n*n);
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                idnmatrix[i*n+j]=1;
            }
            else{
                idnmatrix[i*n+j]=0;
            }
        }
    }
    idn.size=n;
    idn.Mtx=idnmatrix;
    return idn;
}

void SubtractMultipleRows(struct StrMtx inputmatrix,double Coefficient, int lineNumber1, int lineNumber2){
    int i=0;
    for(i=0;i<inputmatrix.width; i++){
        inputmatrix.Mtx[lineNumber1*inputmatrix.width+i]-=Coefficient*inputmatrix.Mtx[lineNumber2*inputmatrix.width+i];
    }
}

void SwapRows(struct StrMtx inputmatrix, int lineNumber1, int lineNumber2){
    int i=0;
    double temp1=0.0;
    double temp2=0.0;
    for(i=0;i<inputmatrix.width; i++){
        temp1=inputmatrix.Mtx[lineNumber1*inputmatrix.width+i];
        temp2=inputmatrix.Mtx[lineNumber2*inputmatrix.width+i];
        inputmatrix.Mtx[lineNumber1*inputmatrix.width+i]=temp2;
        inputmatrix.Mtx[lineNumber2*inputmatrix.width+i]=temp1;
    }

}

void SwapCols(struct StrMtx inputmatrix, int colNumber1, int colNumber2){
    int i=0;
    if(colNumber1==colNumber2){
        return;
    }
    double temp1=0.0;
    double temp2=0.0;
    for(i=0;i<inputmatrix.width; i++){
        temp1=inputmatrix.Mtx[i*inputmatrix.width+colNumber1];
        temp2=inputmatrix.Mtx[i*inputmatrix.width+colNumber2];
        inputmatrix.Mtx[i*inputmatrix.width+colNumber2]=temp1;
        inputmatrix.Mtx[i*inputmatrix.width+colNumber1]=temp2;
    }

}

struct Element GetMaxInRow(struct StrMtx inputmatrix, int lineNumber, int startpoint){
    int i;
    struct Element result;
    result.value=inputmatrix.Mtx[lineNumber*inputmatrix.size+startpoint];
    result.RowPlace=lineNumber;
    result.ColPlace=startpoint;

    for(i=startpoint; i<inputmatrix.size; i++){
        if(fabs(inputmatrix.Mtx[inputmatrix.size*lineNumber+i])> fabs(result.value)){
            result.value=inputmatrix.Mtx[inputmatrix.size*lineNumber+i];
            result.ColPlace=i;
        }
    }
    return result;
}

struct Element GetMaxInSubMtx(struct StrMtx inputmatrix, int lineNumber, int colNumber){
    int i;
    struct Element result;
    struct Element temp;
    result.value=inputmatrix.Mtx[lineNumber*inputmatrix.size+colNumber];
    result.RowPlace=lineNumber;
    result.ColPlace=colNumber;
    temp=result;
    for(i=lineNumber; i<inputmatrix.size; i++){
        temp=GetMaxInRow(inputmatrix, i, colNumber);
        if(fabs(temp.value)>fabs(result.value)){
            result=temp;
        }
    }
    return result;
}

void MultiplyRow(struct StrMtx inputmatrix, int lineNumber, double Coefficient){
    int i;
    for(i=0;i<inputmatrix.width;i++){
        inputmatrix.Mtx[lineNumber*inputmatrix.width+i]=inputmatrix.Mtx[lineNumber*inputmatrix.width+i]*Coefficient;
    }
}

void Norm(struct StrElimMtx matrix){
    int j;
    for(j=0;j<matrix.Left.size;j++){
        double MaxElementValue;
        MaxElementValue=GetMaxInRow(matrix.Left, j, 0).value;
        if(MaxElementValue<threshold){
            printf("A matrix szingularis/rosszul kondicionalt \n");
            exit(0);
        }
        MultiplyRow(matrix.Left, j, 1./(MaxElementValue));
        MultiplyRow(matrix.Right, j, 1./(MaxElementValue));
    }
}

void Pivot(struct StrElimMtx matrix, int rownumber, int* ColPermutation){
    int n=matrix.Left.size;
    if(matrix.Left.Mtx[rownumber*n+rownumber]<threshold){
        struct Element MaxValue=GetMaxInSubMtx(matrix.Left, rownumber, rownumber);
        if(fabs(MaxValue.value)>=threshold){
            SwapRows(matrix.Left, rownumber, MaxValue.RowPlace);
            SwapRows(matrix.Right, rownumber, MaxValue.RowPlace);
            SwapCols(matrix.Left, rownumber, MaxValue.ColPlace);
            SwapCols(matrix.Right, rownumber, MaxValue.ColPlace);
            ColPermutation[rownumber]=MaxValue.ColPlace;
        }

        else{
            printf("A matrix szingularis \n");
            exit(0);
        }
    }


}

void Gauss(struct StrElimMtx matrix){
    int n=matrix.Left.size;
    int k;
    int* permutation;
    permutation=(int *) malloc(n*sizeof(int));
    int d;
    for(d=0;d<n;d++){
        permutation[d]=d;
    }
    Norm(matrix);
    double diag;
    for(k=0;k<n;k++){
        if(fabs(matrix.Left.Mtx[k*n+k])<threshold){
            Pivot(matrix, k, permutation);
        }
        diag=matrix.Left.Mtx[k*n+k];
        int i;
        struct Element Maximum=GetMaxInRow(matrix.Left, k,0);

        if(fabs(Maximum.value)<threshold){
            printf("A matrix szingularis \n");
            exit(0);
        }
        for(i=k+1;i<n;i++){
            double ratio=matrix.Left.Mtx[i*n+k]/diag;
            SubtractMultipleRows(matrix.Left,ratio, i, k);
            SubtractMultipleRows(matrix.Right,ratio, i, k);
        }
        MultiplyRow(matrix.Left, k, 1./diag);
        MultiplyRow(matrix.Right, k, 1./diag);
    }

    for(k=n-1;k>=0;k--){
        int i;
        for(i=k-1;i>=0;i--){
            double ratio=matrix.Left.Mtx[i*n+k];
            SubtractMultipleRows(matrix.Left, ratio, i, k);
            SubtractMultipleRows(matrix.Right, ratio, i, k);
        }
    }
    for(d=0;d<n;d++){
        if(permutation[d]!=d){
            SwapCols(matrix.Left, permutation[d], d);
            SwapRows(matrix.Left, permutation[d], d);
            SwapRows(matrix.Right, permutation[d], d);
        }
    }
}

int main(int argc , char *argv[]){
    if(argc<4){
        printf("%d", argc);
        Usage();
    }
    char* argumentstring1;
    char* argumentstring2;
    if(argc==6){
        argumentstring1=argv[4];
        argumentstring2=argv[5];
    }
    int vars=atoi(argv[1]);
    int deg=atoi(argv[2]);
    char* src=argv[3];
    struct StrMtx argmatrix;
    struct StrMtx error;
    struct IntStrMtx Coeffmatrix=GenerateCoeffs(vars, deg);
    argmatrix=ReadMatrix(src, vars+2);
    error=GenerateErrorVector(argmatrix);
    struct StrMtx terv;
    terv=GenerateTervmtx(argmatrix, Coeffmatrix);
    struct StrMtx a=TransposeMatrix(terv);
    struct StrMtx X=Szorzas(a, terv);
    X.size=X.width;
    struct StrMtx jobbvec=GenerateVector(argmatrix);
    struct StrElimMtx elim;
    elim.Left=X;
    elim.Right=Szorzas(a, jobbvec);
    Gauss(elim);
    FILE* coeffstream;
    if(argc==6){
        coeffstream=fopen(argumentstring1, "w");
    }
    else{
        coeffstream=fopen("egyutthatok_big_harmadfok.txt", "w");
    }
    int i;
    for(i=1;i<vars+1;i++){
        fprintf(coeffstream, "X%i ", i);
        printf("X%i ", i);
    }
    fprintf(coeffstream, "\n");
    printf("\n");
    int j,k;
    for(j=0; j<Coeffmatrix.height; j++){
        for(k=0;k<vars;k++){
            fprintf(coeffstream, "%i  ", Coeffmatrix.Mtx[j*Coeffmatrix.width+k]);
            printf("%i  ", Coeffmatrix.Mtx[j*Coeffmatrix.width+k]);
        }
        fprintf(coeffstream, ":   %lf\n", elim.Right.Mtx[j]);
        printf(":   %lf\n", elim.Right.Mtx[j]);
    }
    fclose(coeffstream);
    struct StrMtx estimate;
    struct StrMtx compareMtx;
    double* compmtx;
    compmtx=(double *)malloc(2*argmatrix.height*sizeof(double));
    compareMtx.height=argmatrix.height;
    compareMtx.width=2;
    estimate=Szorzas(terv, elim.Right);
    double lsq=0.0;
    for(k=0; k<estimate.height; k++){
        compmtx[2*k]=argmatrix.Mtx[k*argmatrix.width+vars];
        compmtx[2*k+1]=estimate.Mtx[k]*error.Mtx[k];
        lsq=lsq+(compmtx[2*k]-compmtx[2*k+1])*(compmtx[2*k]-compmtx[2*k+1])/(error.Mtx[k]*error.Mtx[k]);
    }
    compareMtx.Mtx=compmtx;
    if(argc==6){
        PrintMatrix(compareMtx, argumentstring2);
    }
    else{
        PrintMatrix(compareMtx, "becsult_ertek_big_harmadfok.txt");
    }
    printf("LSQ: %f\n", lsq );
    return 0;
    }
