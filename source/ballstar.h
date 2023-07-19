#ifndef __BALL_STAR_H
#define __BALL_STAR_H
#include "nodeballstar.h"
#include <fstream>
#include <sstream>
#include <string>
#include <queue>

#include "pca.h"

class nearestNeighbor{
public:
    Eigen::MatrixXd point;
    double distance;
    nearestNeighbor(Eigen::MatrixXd _point, double _distance): point(_point), distance(_distance){}

    bool operator<(const nearestNeighbor & _point2) const{
        return distance < _point2.distance;  
    }
};


class BallStar{
public:
    NodeBallStar* root;
    Eigen::MatrixXd data; // conjunto de datos
    BallStar(const string& filename);


    void buildtree();
    NodeBallStar* buildSubtree(Eigen::MatrixXd load_data);
    void split( Eigen::MatrixXd load_data,Eigen::MatrixXd &dataL, Eigen::MatrixXd &dataR );
    


    double bestTC( Eigen::MatrixXd transformada );
    void dividir(Eigen::MatrixXd &dataL, Eigen::MatrixXd &dataR, Eigen::MatrixXd load_data, Eigen::MatrixXd transformada, double bestTc );
    // void CalculateHyperplane( Eigen::MatrixXd load_data, Eigen::MatrixXd transformada, Eigen::MatrixXd& dataL, Eigen::MatrixXd& dataR );




    int sections(int N);
    Eigen::VectorXd getCenter(Eigen::MatrixXd load_data);
    double getRadious(Eigen::MatrixXd load_data, Eigen::MatrixXd center);

    double calculateDistance(const Eigen::VectorXd point1,const Eigen::VectorXd point2);
    Eigen::MatrixXd setdata(const std::string& filename);



    // void searchKNN(vector<Song<N>>& _kNN, Song<N> & _keyPoint, int _k);
    void searchKNN(Eigen::MatrixXd& _KNN, Eigen::MatrixXd _keyPoint, int _k);

    void KNN(priority_queue<nearestNeighbor>& _neighbors, NodeBallStar* _currentNode, Eigen::MatrixXd& _keyPoint, int _k);
    double euclideanDistance(Eigen::MatrixXd& _p1, Eigen::MatrixXd& _p2);

};


#endif
