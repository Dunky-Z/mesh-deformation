#include "stdafx.h"

#include "LaplacianDeform.h"

using namespace std;
using namespace Eigen;
using namespace surface_mesh;

LaplaceDeformation::LaplaceDeformation()
{
}

//初始化：设置固定锚点和移动锚点已经移动锚点移动后的坐标
void LaplaceDeformation::InitializeMesh(const Surface_mesh& mesh)
{
	fix_anchor_idx.clear();
	move_anchor_idx.clear();

	/*fix_anchor_idx.push_back(6);
	fix_anchor_idx.push_back(7);

	move_anchor_coord.resize(mesh.vertices_size());

	move_anchor_idx.push_back(0);
	move_anchor_idx.push_back(1);

	move_anchor_coord[0][0] = -0.7;
	move_anchor_coord[0][1] = -0.7;
	move_anchor_coord[0][2] = 0.7;

	move_anchor_coord[1][0] = 0.7;
	move_anchor_coord[1][1] = -0.7;
	move_anchor_coord[1][2] = 0.7;*/
/////////////////////////////////////////////////////
	/*sphere*/
	/*fix_anchor_idx.push_back(2569);
	fix_anchor_idx.push_back(2570);
	fix_anchor_idx.push_back(2568);

	move_anchor_coord.resize(mesh.vertices_size());

	move_anchor_idx.push_back(1778);
	move_anchor_idx.push_back(4050);

	move_anchor_coord[1778][0] = -0.37;
	move_anchor_coord[1778][1] = 0.359;
	move_anchor_coord[1778][2] = -0.6;

	move_anchor_coord[4050][0] = -0.4;
	move_anchor_coord[4050][1] = 0.33;
	move_anchor_coord[4050][2] = -0.58;	*/
	////////////////////////////////////////////////////////
	//fix_anchor_idx.push_back(671);
	//fix_anchor_idx.push_back(567);
	//fix_anchor_idx.push_back(793);
	//fix_anchor_idx.push_back(646);
	//fix_anchor_idx.push_back(841);
	//fix_anchor_idx.push_back(709);
	//fix_anchor_idx.push_back(396);
	//fix_anchor_idx.push_back(868);
	//fix_anchor_idx.push_back(873);
	//fix_anchor_idx.push_back(967);
	//fix_anchor_idx.push_back(981);
	//fix_anchor_idx.push_back(976);

	//move_anchor_coord.resize(mesh.vertices_size());

	//move_anchor_idx.push_back(1079);
	//move_anchor_idx.push_back(1078);
	//move_anchor_idx.push_back(1085);

	//move_anchor_coord[1079][0] = -484;
	//move_anchor_coord[1079][1] = 241;
	//move_anchor_coord[1079][2] = -4.2;

	//move_anchor_coord[1078][0] = -481.47;
	//move_anchor_coord[1078][1] = 233.41;
	//move_anchor_coord[1078][2] = 9.1;

	//move_anchor_coord[1085][0] = -478.360352;
	//move_anchor_coord[1085][1] = 230.731155;
	//move_anchor_coord[1085][2] = -4.088119;
	/////////////////////////////////////////////////////////////////////
	fix_anchor_idx.push_back(333);
	fix_anchor_idx.push_back(521);
	fix_anchor_idx.push_back(513);
	fix_anchor_idx.push_back(502);
	fix_anchor_idx.push_back(310);
	fix_anchor_idx.push_back(487);

	move_anchor_coord.resize(mesh.vertices_size());

	move_anchor_idx.push_back(498);
	move_anchor_idx.push_back(235);
	move_anchor_idx.push_back(33);
	move_anchor_idx.push_back(258);
	move_anchor_idx.push_back(330);
	move_anchor_idx.push_back(223);
	move_anchor_idx.push_back(395);


	move_anchor_coord[498][0] = 1;
	move_anchor_coord[498][1] = 113.56;
	move_anchor_coord[498][2] = -179.98;
	
	move_anchor_coord[235][0] = -4.02;
	move_anchor_coord[235][1] = 123.56;
	move_anchor_coord[235][2] = -164.84;

	move_anchor_coord[33][0] = -3.71;
	move_anchor_coord[33][1] = 133.151;
	move_anchor_coord[33][2] = -137.28;
	
	move_anchor_coord[258][0] = 0;
	move_anchor_coord[258][1] = 113.151;
	move_anchor_coord[258][2] = -154.28;

	move_anchor_coord[330][0] = 0;
	move_anchor_coord[330][1] = 103.151;
	move_anchor_coord[330][2] = -137.4664;

	move_anchor_coord[223][0] = 3.021;
	move_anchor_coord[223][1] = 110.151;
	move_anchor_coord[223][2] = -154.84;

	move_anchor_coord[395][0] = 3.71;
	move_anchor_coord[395][1] = 120.151;
	move_anchor_coord[395][2] = -137.28;
}

//主函数
void LaplaceDeformation::AllpyLaplaceDeformation(char ** argv)
{
	//获取数据
	Surface_mesh mesh;
	mesh.read(argv[1]);
	cout << "获取数据完成" << endl;

	if (mesh.n_vertices() == 0)
		return;
	InitializeMesh(mesh);
	cout << "初始化完成" << endl;

	BuildAdjacentMatrix(mesh);
	cout << "邻接矩阵" << endl;

	//for (int i=0;i<mesh.n_vertices();i++)
	//{
	//	for (int j = 0; j < mesh.n_vertices(); j++)
	//	{
	//		cout << AdjacentVertices.coeffRef(i,j);
	//	}
	//	cout << endl;
	//}

	//cout << VerticesDegree << endl;

	BuildATtimesAMatrix(mesh);
	cout << "ATA" << endl;

	BuildATtimesbMatrix(mesh);
	cout << "ATb" << endl;

	//通过cholesky分解，解线性方程组 Ax = b，即 ATA * x = LsTb
	//Cholesky分解是把一个对称正定(symmetric, positive definite)的矩阵表示成一个下三角矩阵 L 和其转置 LT 的乘积的分解
	//需要头文件 #include <Eigen/Cholesky>

	//X = ATA.ldlt().solve(B) //X: n*3, B:N*3
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(ATA);
	vx_new = chol.solve(ATbx);
	vy_new = chol.solve(ATby);
	vz_new = chol.solve(ATbz);

	//for (int i = 0; i < mesh.n_vertices(); i++)
	//{
	//	cout << vx_new.coeffRef(i) << endl;
	//}
	//cout << vx_new << endl;
	//cout << vy_new << endl;
	//cout << vz_new << endl;
	SetNewcord(mesh);
	cout << "更新" << endl;

	mesh.write(argv[2]);
}

//建立邻接矩阵和度数矩阵
void LaplaceDeformation::BuildAdjacentMatrix(const Surface_mesh & mesh)
{
	int points_num = mesh.vertices_size();
	AdjacentVertices.resize(points_num, points_num);
	VerticesDegree.resize(points_num, points_num);
	AdjacentVertices.reserve(points_num / 2);
	VerticesDegree.reserve(VectorXi::Constant(points_num, 1));

	for (const auto &v_it : mesh.vertices())
	{
		double VerticesCount = 0.0;

		Surface_mesh::Vertex_around_vertex_circulator vv_it(&mesh, v_it), vv_end(&mesh, v_it);
		do
		{
			AdjacentVertices.insert(v_it.idx(), (*vv_it).idx()) = 1.0;
			VerticesCount++;
		} while (++vv_it != vv_end);
		VerticesDegree.insert(v_it.idx(), v_it.idx()) = VerticesCount;
	}
}

void LaplaceDeformation::BuildATtimesAMatrix(const Surface_mesh & mesh)
{
	const int n_fix_anchors = fix_anchor_idx.size(), n_move_anchors = move_anchor_idx.size(), points_num = mesh.vertices_size();

	L = VerticesDegree - AdjacentVertices;

	A = VerticesDegree - AdjacentVertices;

	//将A矩阵扩展而保持原有数据不变
	A.conservativeResize(points_num + n_fix_anchors + n_move_anchors, points_num);
	
	for (auto i = 0; i < n_fix_anchors; i++)
	{
		for (auto j = 0; j < points_num; j++)
		{
			if (j == fix_anchor_idx[i])
				A.coeffRef(points_num + i, j) = 1;
			else
				A.coeffRef(points_num + i, j) = 0;
		}
	}

	// 移动锚点
	for (auto i = 0; i < move_anchor_idx.size(); i++)
	{
		for (auto j = 0; j < points_num; j++)
		{
			if (j == move_anchor_idx[i])
				A.coeffRef(points_num + n_fix_anchors + i, j) = 1;
			else
				A.coeffRef(points_num + n_fix_anchors + i, j) = 0;
		}
	}
	//printf("Laplace 矩阵计算完成\n");
	//cout << A << endl;
	AT = A.transpose();
	//cout << AT << endl;
	ATA = AT*(A);
	//cout << ATA << endl;
}

void LaplaceDeformation::BuildATtimesbMatrix(const Surface_mesh & mesh)
{

	const int n_fix_anchors = fix_anchor_idx.size(), n_move_anchors = move_anchor_idx.size(), points_num = mesh.vertices_size();
	SparseVector<double> vx(points_num), vy(points_num), vz(points_num);
	int i = 0;
	vector<Point> Vertice;

	for (const auto &v_it : mesh.vertices())
	{
		Point t = mesh.position(v_it);
		Vertice.push_back(t);
		vx.coeffRef(i) = t[0];
		vy.coeffRef(i) = t[1];
		vz.coeffRef(i) = t[2];
		i++;
	}
	// 根据laplace矩阵计算出所有点的的laplace坐标
	bx = L*(vx);
	by = L*(vy);
	bz = L*(vz);

	bx.conservativeResize(points_num + n_fix_anchors + n_move_anchors);
	by.conservativeResize(points_num + n_fix_anchors + n_move_anchors);
	bz.conservativeResize(points_num + n_fix_anchors + n_move_anchors);

	// 用形变前坐标对固定锚点坐标进行赋值
	for (auto i = 0; i < n_fix_anchors; i++)
	{

		bx.insert(i + points_num) = Vertice[fix_anchor_idx[i]][0];
		by.insert(i + points_num) = Vertice[fix_anchor_idx[i]][1];
		bz.insert(i + points_num) = Vertice[fix_anchor_idx[i]][2];
	}
	// 用形变后坐标对移动锚点坐标进行赋值
	for (auto i = 0; i < n_move_anchors; i++)
	{
		bx.insert(i + points_num + n_fix_anchors) = move_anchor_coord[move_anchor_idx[i]][0];
		by.insert(i + points_num + n_fix_anchors) = move_anchor_coord[move_anchor_idx[i]][1];
		bz.insert(i + points_num + n_fix_anchors) = move_anchor_coord[move_anchor_idx[i]][2];
	}
	// 计算三个轴上的 ATb 向量
	ATbx = AT.operator*(bx);
	ATby = AT.operator*(by);
	ATbz = AT.operator*(bz);

}

//更新坐标
void LaplaceDeformation::SetNewcord(Surface_mesh & mesh)
{
	Surface_mesh::Vertex_property<Point> pts = mesh.get_vertex_property<Point>("v:point");
	int i = 0;
	for (const auto &v_it : mesh.vertices())
	{
		if (find(move_anchor_idx.begin(), move_anchor_idx.end(), i) != move_anchor_idx.end())
		{
			pts[v_it][0] = move_anchor_coord[i][0];
			pts[v_it][1] = move_anchor_coord[i][1];
			pts[v_it][2] = move_anchor_coord[i][2];
		}
		else if (find(fix_anchor_idx.begin(), fix_anchor_idx.end(), i) == fix_anchor_idx.end())
		{
			pts[v_it][0] = vx_new.coeffRef(i);
			pts[v_it][1] = vy_new.coeffRef(i);
			pts[v_it][2] = vz_new.coeffRef(i);
		}
		//cout << vx_new.coeffRef(i) << endl;
		i++;

		/*	cout << pts[v_it][0] <<" ";
			cout << pts[v_it][1] <<" ";
			cout << pts[v_it][2] << endl;*/
	}
}
