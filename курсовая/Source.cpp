#include <pcl/visualization/cloud_viewer.h>
#include <iostream>
#include <pcl/io/io.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree_search.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common.h>
#include <pcl/point_types.h>
#include <pcl/octree/octree_pointcloud_voxelcentroid.h>
#include <iostream>
#include <fstream>
#include <pcl/common/transforms.h>
#include <pcl/common/centroid.h>
#include <pcl/common/pca.h>
#include <pcl/filters/voxel_grid.h>
#include <vtkPlaneSource.h>

Eigen::Vector4f PlanebyNormalandPoint(pcl::PointXYZ& Normal, pcl::PointXYZ& Point)//A(x-x0)+B(y-y0)+C(z-z0)=0
{
	Eigen::Vector4f equation;
	equation[0] = Normal.x;//A
	equation[1] = Normal.y;//B
	equation[2] = Normal.z;//C
	equation[3] = (Normal.x * (-Point.x)) + (Normal.y * (-Point.y)) + (Normal.z * (-Point.z));//D
	return equation;
}
float ComputeEq(Eigen::Vector4f& plane, pcl::PointXYZ& Point)
{
	float eq = (plane[0] * Point.x) + (plane[1] * Point.y) + (plane[2] * Point.z) + plane[3];

	return eq;
}
Eigen::Vector4f NormalPlane(Eigen::Vector4f& plane)
{
	Eigen::Vector4f equation;
	float normmult = 1 / (sqrt(plane[0] * plane[0] + plane[1] * plane[1] + plane[2] * plane[2]));
	equation[0] = plane[0] * normmult;//A
	equation[1] = plane[1] * normmult;//B
	equation[2] = plane[2] * normmult;//C
	equation[3] = 0;//D
	return equation;
}

float PointToPointDistance(pcl::PointXYZ& Point1, pcl::PointXYZ& Point2)
{
	float res = sqrt((Point1.x - Point2.x) * (Point1.x - Point2.x) + (Point1.y - Point2.y) * (Point1.y - Point2.y) + (Point1.z - Point2.z) * (Point1.z - Point2.z));

	return res;
}


vtkSmartPointer<vtkPolyData> createPlane(const pcl::ModelCoefficients& coefficients, float scale[2] = nullptr)
{
	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();

	plane->SetNormal(coefficients.values[0], coefficients.values[1], coefficients.values[2]);
	double norm_sqr = coefficients.values[0] * coefficients.values[0]
		+ coefficients.values[1] * coefficients.values[1]
		+ coefficients.values[2] * coefficients.values[2];


	plane->Push(-coefficients.values[3] / sqrt(norm_sqr));
	plane->SetResolution(200, 200);
	plane->Update();

	double pt1[3], pt2[3], orig[3], center[3];
	plane->GetPoint1(pt1);
	plane->GetPoint2(pt2);
	plane->GetOrigin(orig);
	plane->GetCenter(center);

	double _pt1[3], _pt2[3];
	float scale1 = 3.0;
	float scale2 = 3.0;
	if (scale != nullptr)
	{
		scale1 = scale[0];
		scale2 = scale[1];
	}
	for (int i = 0; i < 3; i++) {
		_pt1[i] = scale1 * (pt1[i] - orig[i]);
		_pt2[i] = scale2 * (pt2[i] - orig[i]);
	}
	for (int i = 0; i < 3; ++i)
	{
		pt1[i] = orig[i] + _pt1[i];
		pt2[i] = orig[i] + _pt2[i];
	}
	plane->SetPoint1(pt1);
	plane->SetPoint2(pt2);



	plane->Update();
	return (plane->GetOutput());
}



vtkSmartPointer<vtkPolyData> /*pcl::visualization::*/createPlane(const pcl::ModelCoefficients& coefficients, double x, double y, double z, float scale[2] = nullptr)
{
	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();


	double norm_sqr = 1.0 / (coefficients.values[0] * coefficients.values[0] +
		coefficients.values[1] * coefficients.values[1] +
		coefficients.values[2] * coefficients.values[2]);

	plane->SetNormal(coefficients.values[0], coefficients.values[1], coefficients.values[2]);
	double t = x * coefficients.values[0] + y * coefficients.values[1] + z * coefficients.values[2] + coefficients.values[3];
	x -= coefficients.values[0] * t * norm_sqr;
	y -= coefficients.values[1] * t * norm_sqr;
	z -= coefficients.values[2] * t * norm_sqr;

	plane->SetCenter(x, y, z);

	{
		double pt1[3], pt2[3], orig[3], center[3];
		plane->GetPoint1(pt1);
		plane->GetPoint2(pt2);
		plane->GetOrigin(orig);
		plane->GetCenter(center);

		float scale1 = 3.0;
		float scale2 = 3.0;
		if (scale != nullptr)
		{
			scale1 = scale[0];
			scale2 = scale[1];
		}
		
		double _pt1[3], _pt2[3];
		for (int i = 0; i < 3; i++) {
			_pt1[i] = scale1 * (pt1[i] - orig[i]);
			_pt2[i] = scale2 * (pt2[i] - orig[i]);
		}
		
		for (int i = 0; i < 3; ++i)
		{
			pt1[i] = orig[i] + _pt1[i];
			pt2[i] = orig[i] + _pt2[i];
		}
		plane->SetPoint1(pt1);
		plane->SetPoint2(pt2);


	}
	plane->Update();

	return (plane->GetOutput());
}

Eigen::Vector4f symLRX(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
{
	float ReflectD = 0.0f;
	float ReflectD2 = 0.0f;
	float OriginD = 0.0f;
	float OriginD2 = 0.0f;
	float minDPCA = VTK_FLOAT_MAX;
	std::vector<pcl::PointXYZ> LSET;
	std::vector<pcl::PointXYZ> RSET;
	std::vector<pcl::PointXYZ> ReflectsetR;
	std::vector<pcl::PointXYZ> ReflectsetL;
	float minD = VTK_FLOAT_MAX;
	Eigen::Vector4f minEq;
	pcl::PointXYZ Normalmin;
	int i = 0;
	for (float x1 = 0.9f; x1 <= 1.1f; x1 += 0.02f)
	{
		for (float y1 = -0.1f; y1 <= 0.1f; y1 += 0.02f)
		{
			for (float z1 = -0.1f; z1 <= 0.1f; z1 += 0.02f)
			{
				Eigen::Vector4f plane = { x1,y1,z1,0 };
				for (size_t i = 0; i < cloud->size(); ++i)
				{
					if (ComputeEq(plane, (*cloud)[i]) > 0.0f)
					{
						RSET.push_back((*cloud)[i]);

					}
					else if (ComputeEq(plane, (*cloud)[i]) < 0.0f)
					{
						LSET.push_back((*cloud)[i]);
					}
				}
				for (size_t i = 0; i < LSET.size(); i++)
				{
					pcl::PointXYZ refpoint;
					refpoint.x = (LSET[i].x * (1 - 2 * plane[0] * plane[0])) - (2 * plane[0] * plane[1]) * LSET[i].y - (2 * plane[0] * plane[2]) * LSET[i].z;
					refpoint.y = (LSET[i].y * (1 - 2 * plane[1] * plane[1])) - (2 * plane[0] * plane[1]) * LSET[i].x - (2 * plane[1] * plane[2]) * LSET[i].z;
					refpoint.z = (LSET[i].z * (1 - 2 * plane[2] * plane[2])) - (2 * plane[0] * plane[2]) * LSET[i].x - (2 * plane[1] * plane[2]) * LSET[i].y;
					ReflectsetR.push_back(refpoint);
				}

				for (size_t i = 0; i < RSET.size(); i++)
				{
					pcl::PointXYZ refpoint;
					refpoint.x = (RSET[i].x * (1 - 2 * plane[0] * plane[0])) - (2 * plane[0] * plane[1]) * RSET[i].y - (2 * plane[0] * plane[2]) * RSET[i].z;
					refpoint.y = (RSET[i].y * (1 - 2 * plane[1] * plane[1])) - (2 * plane[0] * plane[1]) * RSET[i].x - (2 * plane[1] * plane[2]) * RSET[i].z;
					refpoint.z = (RSET[i].z * (1 - 2 * plane[2] * plane[2])) - (2 * plane[0] * plane[2]) * RSET[i].x - (2 * plane[1] * plane[2]) * RSET[i].y;
					ReflectsetL.push_back(refpoint);
				}
				float mind1 = VTK_FLOAT_MAX;
				float mind2 = VTK_FLOAT_MAX;
				float mind3 = VTK_FLOAT_MAX;
				for (size_t i = 0; i < RSET.size(); i++)
				{
					for (size_t j = 0; j < LSET.size(); j++)
					{
						float tempmind = PointToPointDistance(RSET[i], LSET[j]);
						if (tempmind < mind1)
						{
							mind1 = tempmind;
						}
					}
					OriginD += mind1/RSET.size();

					for (size_t j = 0; j < ReflectsetR.size(); j++)
					{
						float tempmind = PointToPointDistance(RSET[i], ReflectsetR[j]);
						if (tempmind < mind2)
						{
							mind2 = tempmind;
						}
					}
					ReflectD += mind2/ ReflectsetR.size();


					
					
					mind1 = VTK_FLOAT_MAX;
					mind2 = VTK_FLOAT_MAX;
					mind3 = VTK_FLOAT_MAX;
					
					
					
					
					//std::cout << abs(R) << endl;
					
				}
				for (size_t i = 0; i < LSET.size(); i++)
				{

					for (size_t j = 0; j < ReflectsetL.size(); j++)
					{
						float tempmind = PointToPointDistance(LSET[i], ReflectsetL[j]);
						if (tempmind < mind3)//поиск минимума
						{
							mind3 = tempmind;
						}
					}

					ReflectD2 += mind3/ ReflectsetL.size();

					for (size_t j = 0; j < RSET.size(); j++)
					{
						float tempmind = PointToPointDistance(RSET[j], LSET[i]);
						if (tempmind < mind1)
						{
							mind1 = tempmind;
						}
					}
					OriginD2 += mind1 / LSET.size();

					mind1 = VTK_FLOAT_MAX;
					mind2 = VTK_FLOAT_MAX;
					mind3 = VTK_FLOAT_MAX;
				}

				
				LSET.clear();
				RSET.clear();
				ReflectsetR.clear();
				ReflectsetL.clear();
				if (std::max(ReflectD,ReflectD2) < minD)
				{
					minD = std::max(ReflectD,ReflectD2);
					minEq = plane;
				}
				if (/*x1>0.98f && x1<1.02f && y1>-0.02f &&y1<0.02f && z1>-0.02f && z1 < 0.02f*/i==665)
				{
					float minDPCA = std::max(ReflectD, ReflectD2);
					std::cout << "dHXPCA=" << minDPCA << std::endl;
				}
				OriginD = 0.0f;
				ReflectD = 0.0f;
				OriginD2 = 0.0f;
				ReflectD2 = 0.0f;
				i++;
			}
		}
	}
	std::cout << "dHX="<<minD << std::endl;
	return minEq;
}
Eigen::Vector4f symLRY(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
{
	float ReflectD = 0.0f;
	float ReflectD2 = 0.0f;
	float OriginD = 0.0f;
	float minDPCA = VTK_FLOAT_MAX;
	std::vector<pcl::PointXYZ> LSET;
	std::vector<pcl::PointXYZ> RSET;
	std::vector<pcl::PointXYZ> ReflectsetR;
	std::vector<pcl::PointXYZ> ReflectsetL;
	float minD = VTK_FLOAT_MAX;
	Eigen::Vector4f minEq;
	pcl::PointXYZ Normalmin;
	int i = 0;
	for (float x1 = -0.1f; x1 <= 0.1f; x1 += 0.02f)
	{
		for (float y1 = 0.9f; y1 <= 1.1f; y1 += 0.02f)
		{
			for (float z1 = -0.1f; z1 <= 0.1f; z1 += 0.02f)
			{
				Eigen::Vector4f plane = { x1,y1,z1,0 };
				for (size_t i = 0; i < cloud->size(); ++i)
				{
					if (ComputeEq(plane, (*cloud)[i]) > 0.0f)
					{
						RSET.push_back((*cloud)[i]);

					}
					else if (ComputeEq(plane, (*cloud)[i]) < 0.0f)
					{
						LSET.push_back((*cloud)[i]);
					}
				}
				for (size_t i = 0; i < LSET.size(); i++)
				{
					pcl::PointXYZ refpoint;
					refpoint.x = (LSET[i].x * (1 - 2 * plane[0] * plane[0])) - (2 * plane[0] * plane[1]) * LSET[i].y - (2 * plane[0] * plane[2]) * LSET[i].z;
					refpoint.y = (LSET[i].y * (1 - 2 * plane[1] * plane[1])) - (2 * plane[0] * plane[1]) * LSET[i].x - (2 * plane[1] * plane[2]) * LSET[i].z;
					refpoint.z = (LSET[i].z * (1 - 2 * plane[2] * plane[2])) - (2 * plane[0] * plane[2]) * LSET[i].x - (2 * plane[1] * plane[2]) * LSET[i].y;
					ReflectsetR.push_back(refpoint);
				}

				for (size_t i = 0; i < RSET.size(); i++)
				{
					pcl::PointXYZ refpoint;
					refpoint.x = (RSET[i].x * (1 - 2 * plane[0] * plane[0])) - (2 * plane[0] * plane[1]) * RSET[i].y - (2 * plane[0] * plane[2]) * RSET[i].z;
					refpoint.y = (RSET[i].y * (1 - 2 * plane[1] * plane[1])) - (2 * plane[0] * plane[1]) * RSET[i].x - (2 * plane[1] * plane[2]) * RSET[i].z;
					refpoint.z = (RSET[i].z * (1 - 2 * plane[2] * plane[2])) - (2 * plane[0] * plane[2]) * RSET[i].x - (2 * plane[1] * plane[2]) * RSET[i].y;
					ReflectsetL.push_back(refpoint);
				}
				float mind1 = VTK_FLOAT_MAX;
				float mind2 = VTK_FLOAT_MAX;
				float mind3 = VTK_FLOAT_MAX;
				for (size_t i = 0; i < RSET.size(); i++)
				{
					for (size_t j = 0; j < LSET.size(); j++)
					{
						float tempmind = PointToPointDistance(RSET[i], LSET[j]);
						if (tempmind < mind1)
						{
							mind1 = tempmind;
						}
					}
					OriginD += mind1 / RSET.size();

					for (size_t j = 0; j < ReflectsetR.size(); j++)
					{
						float tempmind = PointToPointDistance(RSET[i], ReflectsetR[j]);
						if (tempmind < mind2)
						{
							mind2 = tempmind;
						}
					}
					ReflectD += mind2 / ReflectsetR.size();




					mind1 = VTK_FLOAT_MAX;
					mind2 = VTK_FLOAT_MAX;
					mind3 = VTK_FLOAT_MAX;




					//std::cout << abs(R) << endl;

				}
				for (size_t i = 0; i < LSET.size(); i++)
				{

					for (size_t j = 0; j < ReflectsetL.size(); j++)
					{
						float tempmind = PointToPointDistance(LSET[i], ReflectsetL[j]);
						if (tempmind < mind3)//поиск минимума
						{
							mind3 = tempmind;
						}
					}

					ReflectD2 += mind3 / ReflectsetL.size();

					mind1 = VTK_FLOAT_MAX;
					mind2 = VTK_FLOAT_MAX;
					mind3 = VTK_FLOAT_MAX;
				}


				LSET.clear();
				RSET.clear();
				ReflectsetR.clear();
				ReflectsetL.clear();
				if (std::max(ReflectD, ReflectD2) < minD)
				{
					minD = std::max(ReflectD, ReflectD2);
					minEq = plane;
				}
				if (i==665)
				{
					float minDPCA = std::max(ReflectD, ReflectD2);
					std::cout << "dHYPCA=" << minDPCA << std::endl;
				}
				OriginD = 0.0f;
				ReflectD = 0.0f;
				ReflectD2 = 0.0f;
				i++;
			}
		}
	}
	
	std::cout << "dHY=" << minD << std::endl;
	return minEq;
}

Eigen::Vector4f symLRZ(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud)
{
	float ReflectD = 0.0f;
	float ReflectD2 = 0.0f;
	float OriginD = 0.0f;
	float minDPCA = VTK_FLOAT_MAX;
	std::vector<pcl::PointXYZ> LSET;
	std::vector<pcl::PointXYZ> RSET;
	std::vector<pcl::PointXYZ> ReflectsetR;
	std::vector<pcl::PointXYZ> ReflectsetL;
	float minD = VTK_FLOAT_MAX;
	Eigen::Vector4f minEq;
	pcl::PointXYZ Normalmin;
	int i = 0;
	for (float x1 = -0.1f; x1 <= 0.1f; x1 += 0.02f)
	{
		for (float y1 = -0.1f; y1 <= 0.1f; y1 += 0.02f)
		{
			for (float z1 = 0.9f; z1 <= 1.1f; z1 += 0.02f)
			{
				Eigen::Vector4f plane = { x1,y1,z1,0 };
				for (size_t i = 0; i < cloud->size(); ++i)
				{
					if (ComputeEq(plane, (*cloud)[i]) > 0.0f)
					{
						RSET.push_back((*cloud)[i]);

					}
					else if (ComputeEq(plane, (*cloud)[i]) < 0.0f)
					{
						LSET.push_back((*cloud)[i]);
					}
				}
				for (size_t i = 0; i < LSET.size(); i++)
				{
					pcl::PointXYZ refpoint;
					refpoint.x = (LSET[i].x * (1 - 2 * plane[0] * plane[0])) - (2 * plane[0] * plane[1]) * LSET[i].y - (2 * plane[0] * plane[2]) * LSET[i].z;
					refpoint.y = (LSET[i].y * (1 - 2 * plane[1] * plane[1])) - (2 * plane[0] * plane[1]) * LSET[i].x - (2 * plane[1] * plane[2]) * LSET[i].z;
					refpoint.z = (LSET[i].z * (1 - 2 * plane[2] * plane[2])) - (2 * plane[0] * plane[2]) * LSET[i].x - (2 * plane[1] * plane[2]) * LSET[i].y;
					ReflectsetR.push_back(refpoint);
				}

				for (size_t i = 0; i < RSET.size(); i++)
				{
					pcl::PointXYZ refpoint;
					refpoint.x = (RSET[i].x * (1 - 2 * plane[0] * plane[0])) - (2 * plane[0] * plane[1]) * RSET[i].y - (2 * plane[0] * plane[2]) * RSET[i].z;
					refpoint.y = (RSET[i].y * (1 - 2 * plane[1] * plane[1])) - (2 * plane[0] * plane[1]) * RSET[i].x - (2 * plane[1] * plane[2]) * RSET[i].z;
					refpoint.z = (RSET[i].z * (1 - 2 * plane[2] * plane[2])) - (2 * plane[0] * plane[2]) * RSET[i].x - (2 * plane[1] * plane[2]) * RSET[i].y;
					ReflectsetL.push_back(refpoint);
				}
				float mind1 = VTK_FLOAT_MAX;
				float mind2 = VTK_FLOAT_MAX;
				float mind3 = VTK_FLOAT_MAX;
				for (size_t i = 0; i < RSET.size(); i++)
				{
					for (size_t j = 0; j < LSET.size(); j++)
					{
						float tempmind = PointToPointDistance(RSET[i], LSET[j]);
						if (tempmind < mind1)
						{
							mind1 = tempmind;
						}
					}
					OriginD += mind1 / RSET.size();

					for (size_t j = 0; j < ReflectsetR.size(); j++)
					{
						float tempmind = PointToPointDistance(RSET[i], ReflectsetR[j]);
						if (tempmind < mind2)
						{
							mind2 = tempmind;
						}
					}
					ReflectD += mind2 / ReflectsetR.size();




					mind1 = VTK_FLOAT_MAX;
					mind2 = VTK_FLOAT_MAX;
					mind3 = VTK_FLOAT_MAX;




					//std::cout << abs(R) << endl;

				}
				for (size_t i = 0; i < LSET.size(); i++)
				{

					for (size_t j = 0; j < ReflectsetL.size(); j++)
					{
						float tempmind = PointToPointDistance(LSET[i], ReflectsetL[j]);
						if (tempmind < mind3)//поиск минимума
						{
							mind3 = tempmind;
						}
					}

					ReflectD2 += mind3 / ReflectsetL.size();

					
					mind1 = VTK_FLOAT_MAX;
					mind2 = VTK_FLOAT_MAX;
					mind3 = VTK_FLOAT_MAX;
				}


				LSET.clear();
				RSET.clear();
				ReflectsetR.clear();
				ReflectsetL.clear();
				if (std::max(ReflectD, ReflectD2) < minD)
				{
					minD = std::max(ReflectD, ReflectD2);
					minEq = plane;
				}
				if (i==665)
				{

					float minDPCA = std::max(ReflectD, ReflectD2);
					std::cout << "dHZPCA=" << minDPCA << std::endl;
				}
				OriginD = 0.0f;
				ReflectD = 0.0f;
				ReflectD2 = 0.0f;
				i++;
			}
		}
	}
	
	std::cout << "dHZ=" << minD << std::endl;
	return minEq;
}
int main()
{
	cout << "Enter file name with 3d object in ply format" << endl;
	std::string file;
	cin >> file;
	if (file == "")
	{
		cout << "file with 3d object not entered" << endl;
		return -1;
	}
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::io::loadPLYFile(file, *cloud);

	vtkObject::GlobalWarningDisplayOff();
	
	
	/*pcl::VoxelGrid <pcl::PointXYZ> avg;
	avg.setInputCloud(cloud);
	avg.setLeafSize(0.2f, 0.2f, 0.2f);
	avg.filter(*cloud);*/

	auto start = std::chrono::high_resolution_clock::now();
	 //вычисление основных компонент
	Eigen::Vector4f pcaCentroid;
	pcl::compute3DCentroid(*cloud, pcaCentroid);
	Eigen::Matrix3f covariance;
	computeCovarianceMatrixNormalized(*cloud, pcaCentroid, covariance);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
	Eigen::Matrix3f eigenVectorsPCA = eigen_solver.eigenvectors();
	eigenVectorsPCA.col(2) = eigenVectorsPCA.col(0).cross(eigenVectorsPCA.col(1));  
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPCAprojection (new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PCA<pcl::PointXYZ> pca;
	pca.setInputCloud(cloud);
	pca.project(*cloud, *cloudPCAprojection);
	

	 //Преобразование исходного облака в начало координат, где главные компоненты соответствуют осям.
	Eigen::Matrix4f projectionTransform(Eigen::Matrix4f::Identity());
	projectionTransform.block<3, 3>(0, 0) = eigenVectorsPCA.transpose();
	projectionTransform.block<3, 1>(0, 3) = -1.f * (projectionTransform.block<3, 3>(0, 0) * pcaCentroid.head<3>());
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPointsProjected(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::transformPointCloud(*cloud, *cloudPointsProjected, projectionTransform);

	
	//Получение минимальных и максимальных точек преобразованного облака
	pcl::PointXYZ minPoint, maxPoint, minPoint1, maxPoint1;
	pcl::getMinMax3D(*cloudPointsProjected, minPoint, maxPoint);
	const Eigen::Vector3f meanDiagonal = 0.5f * (maxPoint.getVector3fMap() + minPoint.getVector3fMap());

	 //окончательное преобразование
	const Eigen::Quaternionf bboxQuaternion(eigenVectorsPCA); //Квартернионы способ вращения
	const Eigen::Vector3f bboxTransform = eigenVectorsPCA * meanDiagonal + pcaCentroid.head<3>();

	
	pcl::visualization::PCLVisualizer visu;
	pcl::visualization::PCLVisualizer visuY;
    pcl::visualization::PCLVisualizer visuZ;
	visu.setBackgroundColor(255,255,255);
	visuY.setBackgroundColor(255, 255, 255);
	visuZ.setBackgroundColor(255, 255, 255);
	visu.setWindowName("Object after normalization");
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color(cloudPointsProjected, 0, 0, 0);
	visu.addPointCloud<pcl::PointXYZ>(cloudPointsProjected,single_color);
	visuY.addPointCloud<pcl::PointXYZ>(cloudPointsProjected, single_color);
	visuZ.addPointCloud<pcl::PointXYZ>(cloudPointsProjected, single_color);
	pcl::VoxelGrid <pcl::PointXYZ> avg;
	avg.setInputCloud(cloudPointsProjected);
	avg.setLeafSize(0.2f, 0.2f, 0.2f);
	avg.filter(*cloudPointsProjected);
	auto end = std::chrono::high_resolution_clock::now();
	auto time = std::chrono::duration<double>(end - start).count();
	pcl::ModelCoefficients plane_coeff;
	Eigen::Vector4f minEq;
	pcl::ModelCoefficients cube_coeff;
	pcl::ModelCoefficients pclX;
	pclX.values.resize(4);   
	pclX.values[0] = 1;
	pclX.values[1] = 0;
	pclX.values[2] = 0;
	pclX.values[3] = 0;
	pcl::ModelCoefficients pclY;
	pclY.values.resize(4);    
	pclY.values[0] = 0;
	pclY.values[1] = 1;
	pclY.values[2] = 0;
	pclY.values[3] = 0;
	pcl::ModelCoefficients pclZ;
	pclZ.values.resize(4);    
	pclZ.values[0] = 0;
	pclZ.values[1] = 0;
	pclZ.values[2] = 1;
	pclZ.values[3] = 0;
	
	pcl::getMinMax3D(*cloudPointsProjected, minPoint1, maxPoint1);
	

	
	
	
	start = std::chrono::high_resolution_clock::now();
	minEq = symLRX(cloudPointsProjected);
	end = std::chrono::high_resolution_clock::now();
	time += std::chrono::duration<double>(end - start).count();
	std::cout << minEq[0] << " " << minEq[1] << " " << minEq[2] << std::endl;
	plane_coeff.values.resize(4);    
	plane_coeff.values[0] = minEq[0];
	plane_coeff.values[1] = minEq[1];
	plane_coeff.values[2] = minEq[2];
	plane_coeff.values[3] = 0;
	
	

	float scale[2] = { 1,2*(maxPoint1.x-minPoint1.x) };
	auto plane = createPlane(plane_coeff, 0, 0, 0, scale);
	auto planex1= createPlane(pclX, 0, 0, 0, scale);
	scale[1] = scale[1];
	scale[0] = -0.5;
	auto planex2 = createPlane(pclX, 0, 0, 0, scale);
	auto plane2= createPlane(plane_coeff, 0, 0, 0, scale);
	
	visu.addModelFromPolyData(plane, "plane_1");
	visu.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1, 0, "plane_1", 0);
	visu.addModelFromPolyData(plane2, "plane_2");
	visu.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1, 0, "plane_2", 0);
	visu.addModelFromPolyData(planex1, "plane_1x");
	visu.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 0, 0, "plane_1x", 0);
	visu.addModelFromPolyData(planex2, "plane_2x");
	visu.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 0, 0, "plane_2x", 0);

	//ось Y
	start = std::chrono::high_resolution_clock::now();
	minEq = symLRY(cloudPointsProjected);
	end = std::chrono::high_resolution_clock::now();
	time += std::chrono::duration<double>(end - start).count();
	std::cout << minEq[0] << " " << minEq[1] << " " << minEq[2] << std::endl;
	plane_coeff.values.resize(4);    
	plane_coeff.values[0] = minEq[0];
	plane_coeff.values[1] = minEq[1];
	plane_coeff.values[2] = minEq[2];
	plane_coeff.values[3] = 0;

	
	scale[0] = 1;
	scale[1] = -(maxPoint1.y-minPoint1.y);
	auto planeY = createPlane(plane_coeff, 0, 0, 0, scale);
	auto planex1Y = createPlane(pclY, 0, 0, 0, scale);
	scale[1] = (maxPoint1.y - minPoint1.y)*1.5;
	auto planex2Y = createPlane(pclY, 0, 0, 0, scale);
	auto plane2Y = createPlane(plane_coeff, 0, 0, 0, scale);
	visuY.addModelFromPolyData(planeY, "plane_1y");
	visuY.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1, 0, "plane_1y", 0);
	visuY.addModelFromPolyData(plane2Y, "plane_2y");
	visuY.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1, 0, "plane_2y", 0);
	visuY.addModelFromPolyData(planex1Y, "plane_1xy");
	visuY.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 0, 0, "plane_1xy", 0);
	visuY.addModelFromPolyData(planex2Y, "plane_2xy");
	visuY.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 0, 0, "plane_2xy", 0);

	//ось Z
	start = std::chrono::high_resolution_clock::now();
	minEq = symLRZ(cloudPointsProjected);
	end = std::chrono::high_resolution_clock::now();
	time += std::chrono::duration<double>(end - start).count();
	std::cout << minEq[0] << " " << minEq[1] << " " << minEq[2] << std::endl;
	std::cout << time;
	plane_coeff.values.resize(4); 
	plane_coeff.values[0] = minEq[0];
	plane_coeff.values[1] = minEq[1];
	plane_coeff.values[2] = minEq[2];
	plane_coeff.values[3] = 0;
	
	scale[0] = 0.75;
	scale[1] = -(maxPoint1.z - minPoint1.z) * 0.15;
	planeY = createPlane(plane_coeff, 0, 0, 0, scale);
	planex1Y = createPlane(pclZ, 0, 0, 0, scale);
	scale[1] = (maxPoint1.z - minPoint1.z)*0.65;
	planex2Y = createPlane(pclZ, 0, 0, 0, scale);
	plane2Y = createPlane(plane_coeff, 0, 0, 0, scale);
	visuZ.addModelFromPolyData(planeY, "plane_1y");
	visuZ.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1, 0, "plane_1y", 0);
	visuZ.addModelFromPolyData(plane2Y, "plane_2y");
	visuZ.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1, 0, "plane_2y", 0);
	visuZ.addModelFromPolyData(planex1Y, "plane_1xy");
	visuZ.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 0, 0, "plane_1xy", 0);
	visuZ.addModelFromPolyData(planex2Y, "plane_2xy");
	visuZ.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 0, 0, "plane_2xy", 0);

	visu.spin();
	visuY.spin();
	visuZ.spin();
	
	return 0;
}