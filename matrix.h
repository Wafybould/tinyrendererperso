#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <vector>
#include <cmath>

class matrix33 {
public:
	std::vector<float> values;
	matrix33() : values(9, 0){}
	matrix33(float m11, float m12, float m13, float m21, float m22, float m23, float m31, float m32, float m33): values(9, 0){
		values[0] = m11;
		values[1] = m12;
		values[2] = m13;
		values[3] = m21;
		values[4] = m22;
		values[5] = m23;
		values[6] = m31;
		values[7] = m32;
		values[8] = m33;
	}
	matrix33(const matrix33& m): values(9, 0){
		values[0] = m.values[0];
		values[1] = m.values[1];
		values[2] = m.values[2];
		values[3] = m.values[3];
		values[4] = m.values[4];
		values[5] = m.values[5];
		values[6] = m.values[6];
		values[7] = m.values[7];
		values[8] = m.values[8];
	}

	float operator[](int i){
		return values[i];
	}

	matrix33& operator*=(float d) {
		for (int i = 0; i < 3; i++) 
		{
			for (int j = 0; j < 3; j++)
			{
				values[i * 3 + j] *= d;
			}
		}
		return *this;
	}

	float determinant() {
		return values[0] * (values[4] * values[8] - values[5] * values[7]) - values[1] * (values[3] * values[8] - values[5] * values[6]) + values[2] * (values[3] * values[7] - values[4] * values[6]);
	}

	float determinant(std::vector<float> v) {
		return v[0] * v[3] - v[1] * v[2];
	}

	matrix33 transpose(matrix33 m) {
		matrix33 res;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				res.values[j * 3 + i] = m.values[i * 3 + j];
			}
		}
		return res;
	}

	matrix33 cofactor(matrix33 m) {
		matrix33 sol;
		std::vector<float> subMat(4);

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				int p = 0;
				for (int k = 0; k < 3; k++)
				{
					if (i == k) {
						continue;
					}
					int q = 0;
					for (int l = 0; l < 3; l++)
					{
						if (j == l) {
							continue;
						}

						subMat[p * 2 + q] = m.values[k * 3 + l];
						q++;
					}
					p++;
				}
				sol.values[i * 3 + j] = pow(-1, i + j) * determinant(subMat);
			}
		}
		return sol;
	}

	matrix33 inverse() {
		float d = 1.0 / determinant();
		matrix33 sol(transpose(cofactor(*this)));
		sol *= d;
		return sol;
	}

	std::vector<float> mult3x1(std::vector<float> m2) {
		std::vector<float> res(3,0);
		res[0] = values[0] * m2[0] + values[1] * m2[1] + values[2] * m2[2];
		res[1] = values[3] * m2[0] + values[4] * m2[1] + values[5] * m2[2];
		res[2] = values[6] * m2[0] + values[7] * m2[1] + values[8] * m2[2];
		return res;
	}

	void print() {
		for (int i = 0; i < 3; i++) {
			std::cout << "(";
			for (int j = 0; j < 3; j++) {
				std::cout << values[i * 3 + j] << "  ";
			}
			std::cout << ")" << std::endl;
		}
	}
};

#endif
