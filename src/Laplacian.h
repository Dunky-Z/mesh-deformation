#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <windows.h>
#include <Eigen/Cholesky>
#include <surface_mesh/Surface_mesh.h>



#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace surface_mesh;
using namespace Eigen;

typedef Eigen::Triplet<double> T;


class LaplaceDeformation
{
public:
	// L: 原始laplace矩阵（方阵），A: 加入了锚点（point_num + anchors, point_num）大小， T：转置
	Eigen::SparseMatrix<double> A, L, AT, ATA, b, ATb, vertice_new;
	std::vector<int>  move_anchor_idx, fix_anchor_idx;
	std::vector<Point> move_anchor_coord;							// 形变后移动锚点位置
	Eigen::SparseMatrix<double> AdjacentVertices, VerticesDegree;	//顶点邻接矩阵;顶点度数矩阵：顶点所相邻顶点的数量

	LaplaceDeformation();
	virtual ~LaplaceDeformation() = default;
	void InitializeMesh();
	void AllpyLaplaceDeformation(char ** argv);
	//void BuildAdjacentMatrix(const Surface_mesh &mesh);
	void Compute_CotMatrix(const Surface_mesh &mesh);
	void BuildATtimesAMatrix(const Surface_mesh & mesh);
	void BuildATtimesbMatrix(const Surface_mesh & mesh);
	void SetNewcord(Surface_mesh & mesh);
};