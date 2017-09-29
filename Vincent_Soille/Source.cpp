// Source.cpp : Defines the entry point for the console application.
//

#include <afxwin.h>  // necessary for MFC to work properly
#include "Header.h"
#include "../../src/blepo.h"
#include <stdlib.h>
#include <iomanip>
#include <stack>
#include <queue>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

using namespace blepo;

// Threshold image
ImgBinary performThreshold(ImgGray img_1, int temp)
{
	ImgBinary img_2;
	img_2.Reset(img_1.Width(), img_1.Height());

	for (int y = 0; y < img_1.Height(); ++y)
	{
		for (int x = 0; x < img_1.Width(); ++x)
		{
			img_2(x, y) = img_1(x, y);
			if (img_1(x, y) > temp)
			{
				img_2(x, y) = 0;
			}
			else
			{
				img_2(x, y) = 1;
			}
		}
	}
	return img_2;
}

// Edge detection
void edgeDetect(ImgInt img_1, ImgBinary *img_2)
{
	(*img_2).Reset(img_1.Width(), img_1.Height());

	for (int y = 0; y < img_1.Height(); ++y)
	{
		for (int x = 0; x < img_1.Width(); ++x)
		{
			(*img_2)(x, y) = 0;
		}
	}

	for (int y = 1; y < img_1.Height() - 1; ++y)
	{
		for (int x = 1; x < img_1.Width() - 1; ++x)
		{
			if (img_1(x, y) != img_1(x + 1, y) || img_1(x, y) != img_1(x - 1, y) || img_1(x, y) != img_1(x, y + 1) || img_1(x, y) != img_1(x, y - 1))
			{
				(*img_2)(x, y) = 1;
			}
		}
	}
}


double calculateKernalHalfWidth(double sigma)
{
	return(floor(2.5*sigma));
}

double calculateKernalWidth(double sigma)
{
	return(floor(2 * (calculateKernalHalfWidth(sigma))) + 1);
}

double* gaussianKernal(double sigma)
{
	double width = calculateKernalWidth(sigma);
	double halfWidth = calculateKernalHalfWidth(sigma);
	double *Gaussian, sum = 0.0;

	Gaussian = new double[(int)width];
	for (int i = 0; i < width; i++)
	{
		Gaussian[i] = exp(-(i - halfWidth)*(i - halfWidth) / (2 * sigma*sigma));
		sum += Gaussian[i];
	}
	for (int i = 0; i < width; i++)
	{
		Gaussian[i] = (Gaussian[i] / sum);
	}

	return Gaussian;
	delete[] Gaussian;
}

double* derivativeKernal(double sigma)
{
	double width = calculateKernalWidth(sigma);
	double halfWidth = calculateKernalHalfWidth(sigma);
	double *Gaussian, sum = 0.0;

	Gaussian = new double[(int)width];
	for (int i = 0; i < width; i++)
	{
		Gaussian[i] = (i - halfWidth)*exp(-(i - halfWidth)*(i - halfWidth) / (2 * sigma*sigma));
		sum += (i*Gaussian[i]);
	}
	for (int i = 0; i < width; i++)
	{
		Gaussian[i] = (Gaussian[i] / sum);
	}

	return Gaussian;
	delete[] Gaussian;
}

ImgFloat convolve(ImgGray img_O, double* temp_1, double* temp_2, double sigma)
{
	ImgFloat img_1;
	img_1.Reset(img_O.Width(), img_O.Height());

	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			img_1(x, y) = img_O(x, y);
		}
	}

	float sum;
	ImgFloat img_2, img_3;
	double width = calculateKernalWidth(sigma);
	double halfWidth = calculateKernalHalfWidth(sigma);

	img_2.Reset(img_O.Width(), img_O.Height());
	img_3.Reset(img_O.Width(), img_O.Height());

	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			img_2(x, y) = 0;
			img_3(x, y) = 0;
		}
	}

	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = halfWidth; x < img_O.Width() - halfWidth; ++x)
		{
			sum = 0;
			for (int z = 0; z < width; z++)
			{
				sum += ((temp_1[z] * img_1(x + halfWidth - z, y)));
			}
		}
	}

	for (int y = halfWidth; y < img_O.Height() - halfWidth; ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			sum = 0;
			for (int z = 0; z < width; z++)
			{
				sum += ((temp_2[z] * img_2(x, y + halfWidth - z)));
			}
			img_3(x, y) = sum;
		}
	}
	return img_3;
}


// Calculate image gradient
void ImgGradient(ImgGray img_O, double sigma, ImgFloat *img_1)
{
	(*img_1).Reset(img_O.Width(), img_O.Height());
	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			(*img_1)(x, y) = 0;
		}
	}

	double *Gaussian = gaussianKernal(sigma);
	double *DGaussian = derivativeKernal(sigma);
	double *flip;
	double width = calculateKernalWidth(sigma);

	flip = new double[(int)width];
	int temp = width - 1;
	for (int z = 0; z < width; z++)
	{
		flip[temp] = DGaussian[z];
		temp--;
	}

	ImgFloat img_2 = convolve(img_O, Gaussian, flip, sigma);
	ImgFloat img_3 = convolve(img_O, flip, Gaussian, sigma);

	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			(*img_1)(x, y) = sqrt((img_3(x, y)*img_3(x, y)) + ((img_2(x, y)*img_2(x, y))));
		}
	}
}

// Calculate chamfer image
void chamferImage(ImgBinary img_1, ImgInt *img_2)
{
	img_2->Reset(img_1.Width(), img_1.Height());
	int temp_1 = (img_1.Width()*img_1.Height()) + 1;

	for (int y = 0; y < img_1.Height(); ++y)
	{
		for (int x = 0; x < img_1.Width(); ++x)
		{
			if (img_1(x, y))
			{
				(*img_2)(x, y) = 0;
			}
			else
			{
				int temp_2 = temp_1;
				if (y > 0)
					temp_2 = blepo_ex::Min(temp_2, (*img_2)(x, y - 1) + 1);
				if (x > 0)
					temp_2 = blepo_ex::Min(temp_2, (*img_2)(x - 1, y) + 1);
				(*img_2)(x, y) = temp_2;
			}
		}
	}

	for (int y = img_1.Height() - 1; y >= 0; y--)
	{
		for (int x = img_1.Width() - 1; x >= 0; x--)
		{
			if (img_1(x, y))
				(*img_2)(x, y) = 0;
			else
			{
				int temp_3 = (*img_2)(x, y);
				if (y < img_1.Height() - 1)
					temp_3 = blepo_ex::Min(temp_3, (*img_2)(x, y + 1) + 1);
				if (x < img_1.Width() - 1)
					temp_3 = blepo_ex::Min(temp_3, (*img_2)(x + 1, y) + 1);
				(*img_2)(x, y) = temp_3;
			}
		}
	}
}

// Compute floodfill
void floodfillImage(ImgInt img_1, int a, int b, int z, ImgInt* img_2)
{
	int temp = img_1(a, b);
	std::stack<Point> temp_1;
	Point temp_2, seedP;
	temp_2.SetPoint(a, b);

	if (temp == z)
		return;

	temp_1.push(Point(a, b));
	(*img_2)(a, b) = z;

	while (temp_1.empty())
	{
		seedP = temp_1.top();
		int x = seedP.x;
		int y = seedP.y;
		temp_1.pop();

		if (img_1(x + 1, y) == temp && (*img_2)(x + 1, y) != z)
		{
			temp_1.push(Point(x + 1, y));
			(*img_2)(x + 1, y) = z;
		}
		if (img_1(x - 1, y) == temp && (*img_2)(x - 1, y) != z)
		{
			temp_1.push(Point(x - 1, y));
			(*img_2)(x - 1, y) = z;
		}
		if (img_1(x + 1, y - 1) == temp && (*img_2)(x + 1, y - 1) != z)
		{
			temp_1.push(Point(x + 1, y - 1));
			(*img_2)(x + 1, y - 1) = z;
		}
		if (img_1(x - 1, y + 1) == temp && (*img_2)(x - 1, y + 1) != z)
		{
			temp_1.push(Point(x - 1, y + 1));
			(*img_2)(x - 1, y + 1) = z;
		}
		if (img_1(x, y + 1) == temp && (*img_2)(x, y + 1) != z)
		{
			temp_1.push(Point(x, y + 1));
			(*img_2)(x, y + 1) = z;
		}
		if (img_1(x, y - 1) == temp && (*img_2)(x, y - 1) != z)
		{
			temp_1.push(Point(x, y - 1));
			(*img_2)(x, y - 1) = z;
		}
		if (img_1(x + 1, y + 1) == temp && (*img_2)(x + 1, y + 1) != z)
		{
			temp_1.push(Point(x + 1, y + 1));
			(*img_2)(x + 1, y + 1) = z;
		}
		if (img_1(x - 1, y - 1) == temp && (*img_2)(x - 1, y - 1) != z)
		{
			temp_1.push(Point(x - 1, y - 1));
			(*img_2)(x - 1, y - 1) = z;
		}
	}
	return;
}

vector<vector<Point>> pixTable(ImgInt img_1)
{
	int size = 0, z = 0;

	for (int y = 0; y < img_1.Height(); ++y)
	{
		for (int x = 0; x < img_1.Width(); ++x)
		{
			if (z < img_1(x, y))
				z = img_1(x, y);
		}
	}

	Point temp;
	temp.SetPoint(0, 0);
	vector<vector<Point>> pTable;

	for (int a = 0; a <= z; a++)
	{
		vector<Point> p;
		for (int y = 0; y < img_1.Height(); ++y)
		{
			for (int x = 0; x < img_1.Width(); ++x)
			{
				if (a == img_1(x, y))
				{
					temp.SetPoint(x, y);
					p.push_back(temp);
				}
			}
		}
		pTable.push_back(p);
	}
	return pTable;
}

//Perform watershed
ImgInt watershedImage(ImgInt img_O)
{
	vector<vector<Point>> pTable;
	pTable = pixTable(img_O);

	ImgInt img_1, img_2;
	img_1.Reset(img_O.Width(), img_O.Height());

	int i, j;

	for (int y = 0; y < img_O.Height(); ++y)
	{
		for (int x = 0; x < img_O.Width(); ++x)
		{
			img_1(x, y) = -1;
		}
	}

	int count = 0;
	std::queue<Point> frontier;
	Point temp;

	for (int a = 0; a < pTable.size(); a++)
	{
		for (int b = 0; b < pTable[a].size(); b++)
		{
			i = pTable[a][b].x;
			j = pTable[a][b].y;

			if (img_1(i + 1, j) != -1 && i<img_1.Width() - 1)
			{
				img_1(i, j) = img_1(i + 1, j);
				temp.SetPoint(i, j);
				frontier.push(temp);
			}
			if (img_1(i - 1, j) != -1 && i>1)
			{
				img_1(i, j) = img_1(i - 1, j);
				temp.SetPoint(i, j);
				frontier.push(temp);
			}
			if (img_1(i, j + 1) != -1 && j<img_1.Height() - 1)
			{
				img_1(i, j) = img_1(i, j + 1);
				temp.SetPoint(i, j);
				frontier.push(temp);
			}
			if (img_1(i, j - 1) != -1 && j>1)
			{
				img_1(i, j) = img_1(i, j - 1);
				temp.SetPoint(i, j);
				frontier.push(temp);
			}
		}

		while (!frontier.empty())
		{
			temp = frontier.front();
			frontier.pop();
			i = temp.x;
			j = temp.y;

			if (img_1(i + 1, j) != -1 && img_O(i + 1, j) <= a && i < img_1.Width() - 1)
			{
				img_1(i + 1, j) = img_1(i, j);
				temp.SetPoint(i, j);
				frontier.push(temp);
			}
			else if (img_1(i - 1, j) != -1 && img_O(i - 1, j) <= a&&i>0)
			{
				img_1(i - 1, j) = img_1(i, j);
				temp.SetPoint(i, j);
				frontier.push(temp);
			}
			else if (img_1(i, j + 1) != -1 && img_O(i, j + 1) <= a&&j<img_1.Height() - 1)
			{
				img_1(i, j + 1) = img_1(i, j);
				temp.SetPoint(i, j);
				frontier.push(temp);
			}
			else if (img_1(i, j - 1) != -1 && img_O(i, j - 1) <= a&&j>0)
			{
				img_1(i, j - 1) = img_1(i, j);
				temp.SetPoint(i, j);
				frontier.push(temp);
			}
		}
	}

	vector<vector<Point>> temp_1;
	temp_1 = pixTable(img_O);

	for (int a = 0; a < temp_1.size(); a++)
	{
		for (int b = 0; b < 10; b++)
		{
			i = temp_1[a][b].x;
			j = temp_1[a][b].y;
			if (img_1(i, j) == -1)
			{
				floodfillImage(img_1, i, j, count, &img_2);
				count++;
			}
		}
	}
	return img_2;
}

int main(int argc, const char* argv[], const char* envp[])
{
	// Initialize MFC and return if failure
	HMODULE hModule = ::GetModuleHandle(NULL);
	if (hModule == NULL || !AfxWinInit(hModule, NULL, ::GetCommandLine(), 0))
	{
		printf("Fatal Error: MFC initialization failed (hModule = %x)\n", (int)hModule);
		return 1;
	}

	int out = 0;

	try
	{
		if (argc > 3)
		{
			cout << "Number of arguments are incorrect" << endl;
			exit(0);
		}
		if (argc < 3)
		{
			cout << "The threshold value for holes.pgm is 77 and cells_small.pgm is 42" << endl;
			exit(0);
		}

		if (argc == 3)
		{
			cout << "The threshold value for holes.pgm is 77 and cells_small.pgm is 42" << endl;
			int out = atoi(argv[2]);

			double sigma = 0.5;

			ImgGray imgFile_Original;

			CString path("../../images/");
			CString imgPath_1 = path + CString(argv[1]);

			Load(imgPath_1, &imgFile_Original);

			// Check for proper file size
			if (imgFile_Original.Width() == 0 || imgFile_Original.Height() == 0)
			{
				printf("Unable to load file 1. Improper file size!\n");
				exit(0);
			}

			// Display the original image
			Figure fig1;
			fig1.SetTitle("Original image");
			fig1.Draw(imgFile_Original);

			ImgFloat imgFile_1;
			ImgGradient(imgFile_Original, sigma, &imgFile_1);
			Figure fig2;
			fig2.SetTitle("Gradient image");
			fig2.Draw(imgFile_1);

			ImgBinary temp;
			temp = performThreshold(imgFile_Original, out);
			Figure fig3;
			fig3.SetTitle("Threshold image");
			fig3.Draw(temp);

			ImgInt imgFile_C;
			chamferImage(temp, &imgFile_C);
			Figure fig4;
			fig4.SetTitle("Chamfer image");
			fig4.Draw(imgFile_C);

			ImgGray imgFile_temp;
			Convert(imgFile_C, &imgFile_temp);

			//watershedImage(imgFile_C); //Unable to resolve this crash issue

			ImgInt imgFile_temp_1;
			ImgBinary imgFile_Edge;
			ImgBinary imgFile_temp_2;

			edgeDetect(imgFile_temp_1, &imgFile_Edge);
			Figure fig5;
			fig5.SetTitle("Edge image");
			fig5.Draw(imgFile_Edge);

			/*Xor(imgFile_Edge, temp, &imgFile_temp_2);
			Figure fig6;
			fig6.SetTitle("Final image");
			fig6.Draw(imgFile_temp_2);*/
			// Loop forever until user presses Ctrl-C in terminal window.
			EventLoop();
		}
	}

	catch (const Exception& e)
	{
		e.Display();    // display exception 
	}
	return 0;
}
