#ifndef __BALL_STAR_H
#define __BALL_STAR_H
#include "nodeballstar.h"
#include <fstream>
#include <sstream>
#include <string>

#include "pca.h"

#define DIMENSIONES 13


class BallStar{
public:
    NodeBallStar* root;
    BallStar();


    void split( const string& filename );
    void split( Eigen::MatrixXd load_data );
    void buildtree(const string& filename);
    void CalculateHyperplane( Eigen::MatrixXd transformada, double N, Eigen::MatrixXd dataL, Eigen::MatrixXd dataR, double radius );







    Eigen::MatrixXd setdata(const std::string& filename);
};


#endif
