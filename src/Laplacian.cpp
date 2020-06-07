
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
	move_anchor_idx.clear();
	/////////////////////////////////////////////////////
	/*sphere*/
	fix_anchor_idx.push_back(2047);

	move_anchor_coord.resize(1);
	move_anchor_idx.push_back(1037);

	move_anchor_coord[0][0] = 0;
	move_anchor_coord[0][1] = 0.490;
	move_anchor_coord[0][2] = -5.99;

	////////////////////////////////////////////////////////
	//fix_anchor_idx.push_back(6243);
	//fix_anchor_idx.push_back(221);
	//fix_anchor_idx.push_back(3906);

	//move_anchor_coord.resize(1);
	//move_anchor_idx.push_back(7276);

	//move_anchor_coord[0][0] = 0.05;
	//move_anchor_coord[0][1] = 0.185;
	//move_anchor_coord[0][2] = -0.0237;


	/////////////////////////////////////////////////////////////////////
	/*long*/
	//fix_anchor_idx.push_back(5044);
	//fix_anchor_idx.push_back(5555);
	//fix_anchor_idx.push_back(5669);
	//fix_anchor_idx.push_back(885);
	//fix_anchor_idx.push_back(5961);
	//fix_anchor_idx.push_back(1066);

	//move_anchor_coord.resize(1);
	//move_anchor_idx.push_back(308);

	//move_anchor_coord[0][0] = 63.5;
	//move_anchor_coord[0][1] = -71;
	//move_anchor_coord[0][2] = 2.1;
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

	Compute_CotMatrix(mesh);
	cout << "余切矩阵计算完成！" << endl;

	/*BuildAdjacentMatrix(mesh);
	cout << "邻接矩阵" << endl;*/

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

	//X = ATA.ldlt().solve(B) //X: n*3, B:N*3
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(ATA);
	vertice_new = chol.solve(ATb);
	cout << "求解完成" << endl;
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


void LaplaceDeformation::Compute_CotMatrix(const Surface_mesh& mesh)
{
	Surface_mesh::Vertex_property<Point> pts = mesh.get_vertex_property<Point>("v:point");
	Surface_mesh::Face_iterator fit = mesh.faces_begin();
	std::vector<T> tripletlist;
	tripletlist.reserve(20);
	const int p_num = mesh.n_vertices();
	L.resize(p_num, p_num);

	do
	{
		Surface_mesh::Vertex_around_face_circulator vf = mesh.vertices(*fit);
		Point p[3];
		double cot[3];
		int id[3];
		for (int i = 0; i < 3; ++i, ++vf)
		{
			p[i] = pts[*vf];
			id[i] = (*vf).idx();
		}

		for (int i = 0; i < 3; ++i)
		{
			int j = (i + 1) % 3, k = (j + 1) % 3;
			cot[i] = dot(p[j] - p[i], p[k] - p[i]) / norm(cross(p[j] - p[i], p[k] - p[i]));

			tripletlist.push_back({ id[j], id[k], -0.5 * cot[i] });
			tripletlist.push_back({ id[k], id[j], -0.5 * cot[i] });
		}
		for (int i = 0; i < 3; ++i)
		{
			tripletlist.push_back({ id[i], id[i], 0.5*(cot[(i + 1) % 3] + cot[(i + 2) % 3]) });
		}

	} while (++fit != mesh.faces_end());




	//int v_id[3];
	//MatrixX3d v(p_num, 3);
	//int k = 0;

	//for (const auto &f_it : mesh.faces())
	//{
	//	k++;
	//	Surface_mesh::Vertex_around_face_circulator vf(&mesh, f_it);
	//	int i = 0;
	//	//获取每个面的三个顶点下标
	//	for (const auto &v_it : vf)
	//	{
	//		v_id[i++] = v_it.idx();
	//	}
	//	int j = 0;
	//	//获取每个顶点的坐标
	//	for (const auto &v_it : vf)
	//	{
	//		v(j, 0) = pts[v_it][0];
	//		v(j, 1) = pts[v_it][1];
	//		v(j, 2) = pts[v_it][2];
	//		j++;
	//	}
	//	double cot1 = (v.row(1) - v.row(0)).dot(v.row(2) - v.row(0)) / ((v.row(1) - v.row(0)).cross(v.row(2) - v.row(0)).norm());
	//	double cot2 = (v.row(2) - v.row(1)).dot(v.row(0) - v.row(1)) / ((v.row(2) - v.row(1)).cross(v.row(0) - v.row(1)).norm());
	//	double cot3 = (v.row(0) - v.row(2)).dot(v.row(1) - v.row(2)) / ((v.row(0) - v.row(2)).cross(v.row(1) - v.row(2)).norm());

	//	tripletlist.push_back(T(v_id[0], v_id[1], -cot3 / 2.0));		tripletlist.push_back(T(v_id[1], v_id[0], -cot3 / 2.0));
	//	tripletlist.push_back(T(v_id[1], v_id[2], -cot1 / 2.0));		tripletlist.push_back(T(v_id[2], v_id[1], -cot1 / 2.0));
	//	tripletlist.push_back(T(v_id[2], v_id[0], -cot2 / 2.0));		tripletlist.push_back(T(v_id[0], v_id[2], -cot2 / 2.0));
	//	tripletlist.push_back(T(v_id[0], v_id[0], (cot2 + cot3) / 2.0));
	//	tripletlist.push_back(T(v_id[1], v_id[1], (cot1 + cot3) / 2.0));
	//	tripletlist.push_back(T(v_id[2], v_id[2], (cot1 + cot2) / 2.0));
	//	/*L.coeffRef(v_id[0], v_id[1]) += -cot3 / 2.0; L.coeffRef(v_id[1], v_id[0]) += -cot3 / 2.0;
	//	L.coeffRef(v_id[1], v_id[2]) += -cot1 / 2.0; L.coeffRef(v_id[2], v_id[1]) += -cot1 / 2.0;
	//	L.coeffRef(v_id[2], v_id[0]) += -cot2 / 2.0; L.coeffRef(v_id[0], v_id[2]) += -cot2 / 2.0;
	//	L.coeffRef(v_id[0], v_id[0]) += cot2 + cot3;
	//	L.coeffRef(v_id[1], v_id[1]) += cot1 + cot3;
	//	L.coeffRef(v_id[2], v_id[2]) += cot1 + cot2;*/

	//}

	L.setFromTriplets(tripletlist.begin(), tripletlist.end());
	//cout << L.topLeftCorner(50, 50) << endl;
	/*for (int _i = 0; _i < p_num; ++_i)
	{
		double sum = 0.0;
		for (int j = 0; j < p_num; ++j)
		{
			sum += L.coeffRef(_i, j);
		}
		L.coeffRef(_i, _i) = -sum;
		cout << "d" << endl;
	}*/
}
////建立邻接矩阵和度数矩阵
//void LaplaceDeformation::BuildAdjacentMatrix(const Surface_mesh & mesh)
//{
//	int points_num = mesh.vertices_size();
//	AdjacentVertices.resize(points_num, points_num);
//	VerticesDegree.resize(points_num, points_num);
//	AdjacentVertices.reserve(points_num / 2);
//	VerticesDegree.reserve(VectorXi::Constant(points_num, 1));
//
//	for (const auto &v_it : mesh.vertices())
//	{
//		double VerticesCount = 0.0;
//		Surface_mesh::Vertex_around_vertex_circulator vv_it(&mesh, v_it), vv_end(&mesh, v_it);
//		do
//		{
//			AdjacentVertices.insert(v_it.idx(), (*vv_it).idx()) = 1.0;
//			VerticesCount++;
//		} while (++vv_it != vv_end);
//		VerticesDegree.insert(v_it.idx(), v_it.idx()) = VerticesCount;
//	}
//}

void LaplaceDeformation::BuildATtimesAMatrix(const Surface_mesh & mesh)
{
	int n_fix_anchors = fix_anchor_idx.size(), n_move_anchors = move_anchor_idx.size(), points_num = mesh.vertices_size();

	//L = VerticesDegree - AdjacentVertices;
	//A = VerticesDegree - AdjacentVertices;

	////将A矩阵扩展而保持原有数据不变
	//A.conservativeResize(points_num + n_fix_anchors + n_move_anchors, points_num);

	//for (auto i = 0; i < n_fix_anchors; i++)
	//{
	//	for (auto j = 0; j < points_num; j++)
	//	{
	//		if (j == fix_anchor_idx[i])
	//			A.coeffRef(points_num + i, j) = 1;
	//		else
	//			A.coeffRef(points_num + i, j) = 0;
	//	}
	//}

	//// 移动锚点
	//for (auto i = 0; i < move_anchor_idx.size(); i++)
	//{
	//	for (auto j = 0; j < points_num; j++)
	//	{
	//		if (j == move_anchor_idx[i])
	//			A.coeffRef(points_num + n_fix_anchors + i, j) = 1;
	//		else
	//			A.coeffRef(points_num + n_fix_anchors + i, j) = 0;
	//	}
	//}
	////printf("Laplace 矩阵计算完成\n");
	////cout << A << endl;
	//AT = A.transpose();
	////cout << AT << endl;
	//ATA = AT * A;
	////cout << ATA << endl;

	/////////////////////////////////////////////////////////////////////////////
	//将A矩阵扩展而保持原有数据不变
	A = L;
	A.conservativeResize(points_num + n_fix_anchors + n_move_anchors, points_num);

	for (auto i = 0; i < n_fix_anchors; i++)
	{
		for (auto j = 0; j < points_num; j++)
		{
			if (j == fix_anchor_idx[i])
				A.coeffRef(points_num + i, j) = 1;
			else
				A.coeffRef(points_num + i, j) = 0;
			//cout << j << endl;
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
			//cout << j << endl;
		}
	}
	//printf("Laplace 矩阵计算完成\n");
	//cout << A << endl;
	AT = A.transpose();
	//cout << AT << endl;
	ATA = AT * A;
	//cout << ATA << endl;
}

void LaplaceDeformation::BuildATtimesbMatrix(const Surface_mesh & mesh)
{
	int n_fix_anchors = fix_anchor_idx.size(), n_move_anchors = move_anchor_idx.size(), points_num = mesh.vertices_size();
	SparseMatrix<double> v(points_num, 3);
	std::vector<T> tripletlist;
	tripletlist.reserve(3);
	int i = 0;
	vector<Point> Vertice;

	for (const auto &v_it : mesh.vertices())
	{
		Point t = mesh.position(v_it);
		Vertice.push_back(t);
		tripletlist.push_back(T(i, 0, t[0]));
		tripletlist.push_back(T(i, 1, t[1]));
		tripletlist.push_back(T(i, 2, t[2]));
		i++;
	}
	v.setFromTriplets(tripletlist.begin(), tripletlist.end());

	// 根据laplace矩阵计算出所有点的的laplace坐标
	b = L * v;
	b.conservativeResize(points_num + n_fix_anchors + n_move_anchors, 3);

	// 用形变前坐标对固定锚点坐标进行赋值
	for (auto i = 0; i < n_fix_anchors; i++)
	{

		b.insert(i + points_num, 0) = Vertice[fix_anchor_idx[i]][0];
		b.insert(i + points_num, 1) = Vertice[fix_anchor_idx[i]][1];
		b.insert(i + points_num, 2) = Vertice[fix_anchor_idx[i]][2];
	}

	// 用形变后坐标对移动锚点坐标进行赋值
	for (auto i = 0; i < n_move_anchors; i++)
	{
		b.insert(i + points_num + n_fix_anchors, 0) = move_anchor_coord[i][0];
		b.insert(i + points_num + n_fix_anchors, 1) = move_anchor_coord[i][1];
		b.insert(i + points_num + n_fix_anchors, 2) = move_anchor_coord[i][2];
	}
	// 计算三个轴上的 ATb 向量
	ATb = AT * b;
}
 
//更新坐标
void LaplaceDeformation::SetNewcord(Surface_mesh & mesh)
{
	Surface_mesh::Vertex_property<Point> pts = mesh.get_vertex_property<Point>("v:point");
	int i = 0;
	for (const auto &v_it : mesh.vertices())
	{
		//在移动锚点了没找到。在固定锚点里没找到则更新
		/*if (find(fix_anchor_idx.begin(), fix_anchor_idx.end(), i) == fix_anchor_idx.end())
		{*/
		pts[v_it][0] = vertice_new.coeffRef(i, 0);
		pts[v_it][1] = vertice_new.coeffRef(i, 1);
		pts[v_it][2] = vertice_new.coeffRef(i, 2);
		//}
		i++;
	}
}
