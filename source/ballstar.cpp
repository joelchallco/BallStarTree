#include "ballstar.h"

using namespace std;

BallStar::BallStar(){
    root = nullptr;
}


Eigen::MatrixXd BallStar::setdata(const string& filename){
    // std::ifstream file(filename);
    // int n=100;
    // Eigen::MatrixXd pca_data_matrix(n, DIMENSIONES);

    std::ifstream file(filename);
    std::vector<std::vector<double>> vec;
    std::string line;
    std::vector<string> name;

    std::getline(file, line);

    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            try {
                double num = std::stod(value);
                row.push_back(num);
            } catch (const exception& e) {
                // Ignore non-numeric values
                name.push_back(value);
            }
        }

        if (!row.empty()) {
            vec.push_back(row);
        }
    }

    Eigen::MatrixXd points(vec.size(), vec[0].size());
    for (int i = 0; i < vec.size(); ++i) {
        for (int j = 0; j < vec[0].size(); ++j) {
            points(i, j) = vec[i][j];
        }
    }

    // cout<<points;
    return points;
}

void BallStar::split(Eigen::MatrixXd load_data){
    
    // Eigen::MatrixXd load_data = setdata(filename);
    // cout<<load_data;

    //Aplica PCA a los puntos para poder obtener los autovectores
	pca_t<double> pca;
	pca.set_input(load_data);
	pca.compute();
    // cout << "Eigen Vectors:		\n" << pca.get_eigen_vectors() << std::endl << std::endl;

    Eigen::MatrixXd eigenVector = pca.get_eigen_vectors();

    // cout<<"EigenVectors"<<eigenVector.rows()<<endl;

    // primer auto vector
    Eigen::VectorXd firstEigenVector = eigenVector.col(0);

    //tranformada
    Eigen::MatrixXd T =  load_data * firstEigenVector ;
    cout<<"T: " <<T<<endl;

    


    // Segundo objetivo,


    // std::vector<Point> leftPoints;
    // std::vector<Point> rightPoints;
    // split(points, center, radius, leftPoints, rightPoints);

    Eigen::MatrixXd dataL;
    Eigen::MatrixXd dataR;
    double radius;
    CalculateHyperplane( T, 35000.0, dataL, dataR, radius );

}


void BallStar::CalculateHyperplane( Eigen::MatrixXd transformada, double N, Eigen::MatrixXd dataL, Eigen::MatrixXd dataR, double radius ){

    double alfa = 0.1;
    double Tmax = transformada.maxCoeff();
    double Tmin = transformada.minCoeff();

    std::cout << "T max: " << Tmax << std::endl;
    std::cout << "T min: " << Tmin << std::endl;


    // 2do objetivo, Minimizar el radio de cada partición de los hijos (R1 y R2)
    double tc;
    double obj2;
    tc = (Tmax + Tmin)/2;
    // cout<<"tc:" <<tc<<endl;

    obj2 =alfa*( (tc-Tmin) / Tmax-Tmin );
    cout<<"obj2: "<<obj2<<endl;

    // 1er objetivo, Maximizar el equilibrio entre el número de puntos en las particiones de los hijos (N1 y N2)
    int tmpIzq=0,tmpDer=0;
    for(int i=0;i<transformada.size();i++){
        if(transformada(i) < tc){
            tmpIzq++;
        }
        else{
            tmpDer++;
        }
    }
    // cout<<"izquierdo: "<<tmpIzq<<endl;
    // cout<<"derecho: "<<tmpDer<<endl;

    double obj1;
    obj1 = abs(tmpIzq - tmpDer)/N;
    cout<<"obj1: "<<obj1<<endl;

    double res = obj1 + obj2;
    cout<<"tc: "<<res<<endl;

    // for(int =0;i<transformada.size();i++){
    //     if()
    // }

}

void BallStar::buildtree(const string& filename){
    // Eigen::MatrixXd load_data = setdata(filename);
    Eigen::MatrixXd load_data = setdata(filename);

    split(load_data);


}



