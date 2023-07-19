#include "ballstar.h"

using namespace std;


BallStar::BallStar(const string& filename){
    root = nullptr;
    data = setdata(filename);
    cout<<"data: "<<data.rows()<<endl<<endl;
}

Eigen::MatrixXd BallStar::setdata(const string& filename){

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

    Eigen::MatrixXd data;
    data.resize(vec.size(), DIMENSIONES);
    for (int i = 0; i < vec.size(); ++i) {
        data(i, 0) = i;
        
        for (int j = 0; j < vec[0].size(); ++j) {
            data(i, j+1) = vec[i][j];

        }
    }
    return data;
}

Eigen::VectorXd BallStar::getCenter(Eigen::MatrixXd load_data){
    int numPoints = load_data.rows(); // Número total de puntos

    Eigen::VectorXd sum = load_data.colwise().sum(); // Suma de las coordenadas
    Eigen::VectorXd center = sum / numPoints; // Cálculo del centro

    // cout << "Centro: " << center.transpose() << std::endl;

    return center.transpose();
}

double BallStar::getRadious(Eigen::MatrixXd load_data, Eigen::MatrixXd center ){
    int numPoints = load_data.rows(); // Número total de puntos
    double radious = 0.0;

    // Calcula la distancia máxima desde el centro a los puntos
    for (int i = 0; i < numPoints; i++) {
        Eigen::VectorXd point = load_data.row(i);
        double distance = calculateDistance(center, point);
        if (distance > radious) {
            radious = distance;
        }
    }

    // std::cout << "Radio: " << radious << std::endl;
    return radious;
}

double BallStar::calculateDistance(const Eigen::VectorXd point1,const Eigen::VectorXd point2) {
    Eigen::VectorXd diff = point1 - point2;
    return diff.norm(); //norma del vector diff
}

double BallStar::euclideanDistance(Eigen::MatrixXd& _p1, Eigen::MatrixXd& _p2){
    
    Eigen::MatrixXd diff = _p1 - _p2;
    // double distance = diff.array().square().sum();
    return diff.lpNorm<2>(); 
    
}


int BallStar::sections(int N){
    int S=0;
    if(N < 40000){ S = 10000; }
    else if(N < 20000){ S = 5000;}
    else if(N < 10000){ S = 2500;}
    else if(N < 5000){ S = 1000;}
    else if(N < 2500){ S = 500;}
    else if(N < 1000){ S = 250;}
    else if(N < 500){ S = 100;}
    else if(N < 250){ S = 50;}
    else if(N < 100){ S = 25;}
    else if(N < 50){ S = 10;}
    // else if(N < 10){ S = 2;}
    return S;
}

void BallStar::dividir(Eigen::MatrixXd &dataL, Eigen::MatrixXd &dataR, Eigen::MatrixXd load_data, Eigen::MatrixXd transformada, double bestTc ){

int indiceR = 0, indiceL = 0;

    dataR.resize(transformada.cols(), load_data.cols()); // Redimensionar dataR para que tenga el mismo número de filas que transformada
    dataL.resize(transformada.cols(), load_data.cols()); // Redimensionar dataL para que tenga el mismo número de filas que transformada

    for (int i = 0; i < transformada.cols(); i++) {
        // cout<<"compara transformada(i) < bestTc : "<< transformada(i) << " - " << bestTc << endl<<endl;
        if (transformada(0,i) < bestTc) {
            dataR.row(indiceR) = load_data.row(i); 
            indiceR++;
        } else {
            dataL.row(indiceL) = load_data.row(i);
            indiceL++;
        }
    }

    dataR.conservativeResize(indiceR, load_data.cols()); // Redimensionar dataR para ajustarse al número real de filas
    dataL.conservativeResize(indiceL, load_data.cols()); // Redimensionar dataL para ajustarse al número real de filas
}

double BallStar::bestTC( Eigen::MatrixXd transformada ){
    double N = transformada.cols();
    int S = sections(N);
    double alfa = 0.1;

    Eigen::Index maxRow, maxCol;
    double tmax = transformada.maxCoeff(&maxRow, &maxCol);
    double tmin = transformada.minCoeff(&maxRow, &maxCol);

    // para saber cuantos puntos estan a la derecha y a la izquierda
    double range = tmax - tmin;
    double tamSection = range/S;
    double besttc; 
    double valorOpt = numeric_limits<double>::max(); 
    double f1,f2;
    double midleSection = tmin + tamSection/2; 
    int N1,N2;

    for(int i=0;i<S;i++){
        N1=0,N2=0;
        f2 = alfa * ((midleSection - tmin) / range );

        // contar N1, N2
        for (int j = 0; j < transformada.cols(); j++){
            if( transformada(0,j) < midleSection){
                N1++;
            }else {
                N2++;
            }  
        }
        // cout << "N1: " << N1 << ", N2: " << N2 << endl;
        f1 = abs(N2 - N1) / ( N*1.0 );
        
        if((f1+f2) < valorOpt){
            valorOpt = f1 + f2;
            besttc = midleSection;
        }
        midleSection += tamSection; // next section
    }
    // cout << "N1: " << N1 << ", N2: " << N2 << endl<<endl;

    return besttc;
}

void BallStar::split(Eigen::MatrixXd load_data, Eigen::MatrixXd &dataL, Eigen::MatrixXd &dataR ){

    //Aplica PCA a los puntos para poder obtener los autovectores
	pca_t<double> pca;
	pca.set_input(load_data);
	pca.compute();

    Eigen::MatrixXd eigenVector = pca.get_eigen_vectors();
    Eigen::VectorXd firstEigenVector = eigenVector.col(0);
    // cout<<"primer autovector"<<firstEigenVector<<endl<<endl;

    Eigen::MatrixXd load_data_centrada = pca.get_centered_matrix();
    load_data_centrada.transposeInPlace();

    //tranformada
    Eigen::MatrixXd transformada = firstEigenVector.transpose() * load_data_centrada;
    
    
    // CalculateHyperplane( load_data, transformada, dataL, dataR );
    double besttc = bestTC(transformada );
    // cout<< "best TC: "<<besttc<<endl<<endl;

    dividir(dataL, dataR,load_data, transformada, besttc);
    

    // cout<<"  Data:"<< load_data.rows()<<endl;
    // cout<<"Data L:"<< dataL.rows()<<endl;
    // cout<<"Data R:"<< dataR.rows()<<endl;
    // cout<<endl;

    // cout<<endl;
}

// void BallStar::CalculateHyperplane( Eigen::MatrixXd load_data, Eigen::MatrixXd transformada, Eigen::MatrixXd& dataL, Eigen::MatrixXd& dataR ){
//     double besttc = bestTC(transformada );
//     // cout<< "best TC: "<<besttc<<endl<<endl;

//     dividir(dataL, dataR,load_data, transformada, besttc);
    

//     // cout<<"  Data:"<< load_data.rows()<<endl;
//     // cout<<"Data L:"<< dataL.rows()<<endl;
//     // cout<<"Data R:"<< dataR.rows()<<endl;
//     // cout<<endl;
// }


void BallStar::buildtree(){

    root = buildSubtree(data);
    
    // cout<<"raiz: "<<endl<<endl;
    // cout<<root->center<<endl;
    // cout<<root->radious<<endl;
    // cout<<root->isleaf<<endl;
}

NodeBallStar* BallStar::buildSubtree(Eigen::MatrixXd load_data){
    if ( load_data.rows() < 138){
        // cout<<"hoja"<<endl;
        // cout<<load_data.rows()<<endl<<endl;
        Eigen::VectorXd center = getCenter(load_data);
        double radious = getRadious(load_data,center);
        NodeBallStar* node = new NodeBallStar(load_data, center.transpose(), radious, true);
        // node->isleaf = true;
        return node;
    }

    // cout<<"---------------------------------------------------------------------------------------------"<<endl;
    Eigen::MatrixXd dataL;
    Eigen::MatrixXd dataR;
    split(load_data, dataL, dataR);

    NodeBallStar* leftNode = buildSubtree(dataL);
    NodeBallStar* rightNode = buildSubtree(dataR);

    Eigen::VectorXd center = getCenter(load_data);
    double radious = getRadious(load_data,center);
    NodeBallStar* node = new NodeBallStar(center.transpose(), radious);
    node->left = leftNode;
    node->right = rightNode;
    // cout<<"center: "<<node->center<<endl;
    // cout<<"radoius: "<<node->radious<<endl;
    return node;
}

void BallStar::KNN(priority_queue<nearestNeighbor>& _neighbors, NodeBallStar* _currentNode, Eigen::MatrixXd& _keyPoint, int _k){
    
    // cout<<"key point:"<<endl<<_keyPoint<<endl<<endl;
    // cout<<"current node:"<<endl<<_currentNode->center<<endl<<endl;
    // cout<<euclideanDistance(_keyPoint, _currentNode->center)<<endl;
    if(euclideanDistance(_keyPoint, _currentNode->center) >= _neighbors.top().distance){
        return; // sin cambios
    }
    
    if (_currentNode->isleaf ){
        for(int i = 0; i < _currentNode->data.rows(); i++){
            Eigen::MatrixXd cancion = _currentNode->data.row(i);
            double tmpDistance = euclideanDistance(_keyPoint, cancion);
            if(tmpDistance < _neighbors.top().distance){
                _neighbors.push(nearestNeighbor(cancion, tmpDistance));
            }
            if(_neighbors.size() > _k){
                _neighbors.pop();
                
            }
        }
    }
    else{
    // left mas cercano, right mas lejano
    NodeBallStar* left;
    NodeBallStar* right;

    if(euclideanDistance(_currentNode->center, _currentNode->left->center) < 
        euclideanDistance(_currentNode->center, _currentNode->right->center) ){
        left = _currentNode->left;
        right = _currentNode->right;
    } else {
        left = _currentNode->right;
        right = _currentNode->left;
    }
    
    KNN(_neighbors, left, _keyPoint, _k);
    KNN(_neighbors, right, _keyPoint, _k);
    }
}



void BallStar::searchKNN(Eigen::MatrixXd& knn, Eigen::MatrixXd _keyPoint, int _k){
    
    priority_queue<nearestNeighbor> priorityQueue;
    Eigen::MatrixXd h;
    h.resize(1,DIMENSIONES);
    // h << 0,0,0,0,0,0,0,0,0,0,0,0,0;
    nearestNeighbor NN( h, numeric_limits<double>::max() );
    priorityQueue.push(NN);
    
    KNN( priorityQueue, root, _keyPoint, _k );
    
    int n = priorityQueue.size();
    knn.resize(n, DIMENSIONES);

    // Copiar los elementos de la priorityQueue a la matriz _KNN
    for (int i = n-1; i >= 0; i--) {
        knn.row(i) = priorityQueue.top().point;
        priorityQueue.pop();
    }
}