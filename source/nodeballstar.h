#ifndef __NODE_BALL_STAR_H
#define __NODE_BALL_STAR_H
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#define DIMENSIONES 14
using namespace std;

class song{
public:
	vector<string> name;
	Eigen::MatrixXd data;
};

class NodeBallStar{
public:
	string name;
    // vector<double> data;
	// dataR.resize(transformada.size(), load_data.cols()); // Redimensionar dataR para que tenga el mismo número de filas que transformada
	Eigen::MatrixXd data;
	// song data;
    Eigen::MatrixXd center;
    double radious;
	bool isleaf;

    NodeBallStar* left;
    NodeBallStar* right;

	NodeBallStar( Eigen::MatrixXd _data, Eigen::MatrixXd _center, double _radious, bool _isleaf ){
		// this->name = _name;
		this->data = _data;
		this->center = _center;
		this->radious = _radious;
		left = nullptr;
		right = nullptr;
		isleaf = _isleaf;
		// center.resize(1,DIMENSIONES);
		// dataR.resize(transformada.size(), load_data.cols()); // Redimensionar dataR para que tenga el mismo número de filas que transformada
		// data.resize(1,DIMENSIONES);
	};

	NodeBallStar( Eigen::MatrixXd _center, double _radious ){
		// this->name = _name;
		this->center = _center;
		this->radious = _radious;
		left = nullptr;
		right = nullptr;
		isleaf = false;
		// dataR.resize(transformada.size(), load_data.cols()); // Redimensionar dataR para que tenga el mismo número de filas que transformada
		// data.resize(1,DIMENSIONES);
	};

};



#endif