#include "Laplacian.h"

using namespace std;
using namespace Eigen;
using namespace surface_mesh;

LaplaceDeformation::LaplaceDeformation()
{
}

//初始化：设置固定锚点和移动锚点已经移动锚点移动后的坐标
void LaplaceDeformation::InitializeMesh()
{
	fix_anchor_idx.clear();
	/*随机生成固定点*/
	int a[4098];
	for (auto i = 0; i <= 4097; ++i) a[i] = i;
	for (auto i = 4097; i >= 1; --i) swap(a[i], a[rand() % i]);
	///////////////////////////////////////////////////
		/*sphere*/
	int id = 0;
	while (id < 150)
	{
		fix_anchor_idx.push_back(a[id]);
		id++;
	}
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
	InitializeMesh();
	cout << "初始化完成" << endl;

	BuildAdjacentMatrix(mesh);
	cout << "邻接矩阵" << endl;


	BuildATtimesAMatrix(mesh);
	cout << "ATA" << endl;

	BuildATtimesbMatrix(mesh);
	cout << "ATb" << endl;

	//通过cholesky分解，解线性方程组 Ax = b，即 ATA * x = ATb
	//Cholesky分解是把一个对称正定(symmetric, positive definite)的矩阵表示成一个下三角矩阵 L 和其转置 LT 的乘积的分解
	//需要头文件 #include <Eigen/Cholesky>

	//X = ATA.ldlt().solve(B) //X: n*3, B:N*3
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(ATA);
	vertice_new = chol.solve(ATb);
	cout << "求解完成" << endl;

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
	VerticesDegree.makeCompressed();
	AdjacentVertices.makeCompressed();
}

void LaplaceDeformation::BuildATtimesAMatrix(const Surface_mesh & mesh)
{
	int n_fix_anchors = fix_anchor_idx.size(), points_num = mesh.vertices_size();

	//L = VerticesDegree - AdjacentVertices;
	A = VerticesDegree - AdjacentVertices;

	//将A矩阵扩展而保持原有数据不变
	A.conservativeResize(points_num + n_fix_anchors, points_num);

	for (auto i = 0; i < n_fix_anchors; i++)
	{
		for (auto j = 0; j < points_num; j++)
		{
			if (j == fix_anchor_idx[i])
				A.insert(points_num + i, j) = 1;
			else
				A.insert(points_num + i, j) = 0;
		}
	}
	A.makeCompressed();
	//printf("Laplace 矩阵计算完成\n");
	//cout << A << endl;
	AT = A.transpose();
	AT.makeCompressed();
	//cout << AT << endl;
	ATA = AT * A;
	ATA.makeCompressed();
	//cout << ATA << endl;
}

void LaplaceDeformation::BuildATtimesbMatrix(const Surface_mesh & mesh)
{
	int n_fix_anchors = fix_anchor_idx.size(), points_num = mesh.vertices_size();
	SparseMatrix<double> v(points_num, points_num);
	int i = 0;
	vector<Point> Vertice;

	for (const auto &v_it : mesh.vertices())
	{
		Point t = mesh.position(v_it);
		Vertice.push_back(t);
		v.coeffRef(i, 0) = t[0];
		v.coeffRef(i, 1) = t[1];
		v.coeffRef(i, 2) = t[2];
		i++;
	}
	// 根据laplace矩阵计算出所有点的的laplace坐标
	b.setZero();

	b.conservativeResize(points_num + n_fix_anchors, 3);

	// 用形变前坐标对固定锚点坐标进行赋值
	for (auto i = 0; i < n_fix_anchors; i++)
	{

		b.insert(i + points_num, 0) = Vertice[fix_anchor_idx[i]][0];
		b.insert(i + points_num, 1) = Vertice[fix_anchor_idx[i]][1];
		b.insert(i + points_num, 2) = Vertice[fix_anchor_idx[i]][2];
	}

	// 计算三个轴上的 ATb 向量
	b.makeCompressed();
	ATb = AT * b;
	ATb.makeCompressed();
}

//更新坐标
void LaplaceDeformation::SetNewcord(Surface_mesh & mesh)
{
	Surface_mesh::Vertex_property<Point> pts = mesh.get_vertex_property<Point>("v:point");
	int i = 0;
	for (const auto &v_it : mesh.vertices())
	{
		pts[v_it][0] = vertice_new.coeffRef(i, 0);
		pts[v_it][1] = vertice_new.coeffRef(i, 1);
		pts[v_it][2] = vertice_new.coeffRef(i, 2);
		i++;
	}
}
