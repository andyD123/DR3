#include "pch.h"

#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/vec_bool_d.h"
#include "../Vectorisation/VecX/vec_double.h"
#include  "../Vectorisation/VecX/alloc_policy.h"
#include  "../Vectorisation/VecX/alloc_policy_imp.h"
#include "../Vectorisation/VecX/target_name_space.h"


TEST(TestCaseAlloc, fillup_empty_last) {
  EXPECT_EQ(1, 1);
  EXPECT_TRUE(true);

  PoolStrat<double> myPool(10);
  int  MAX_EL = 20;
  std::vector<double*> pAlloc;

  int pos = myPool.pos();
  for (int i = 0; i < MAX_EL; ++i)
  {
	  double* p = myPool.alloc();
	  pos = myPool.pos();

	  (*p) = i;
	  pAlloc.push_back(p);
  }

  int szx =myPool.size();
  pos = myPool.pos();

  for (int k = pos; k > 0; --k)
  {
	  double* pback = pAlloc.back();
		  pAlloc.pop_back();
	  myPool.free(pback);

	  pos = myPool.pos();

  }


}


TEST(TestCaseAlloc, fillup_empty_secondlast) {
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	PoolStrat<double> myPool(10);
	int  MAX_EL = 20;
	std::vector<double*> pAlloc;

	int pos = myPool.pos();
	for (int i = 0; i < MAX_EL; ++i)
	{
		double* p = myPool.alloc();
		pos = myPool.pos();

		(*p) = i;
		pAlloc.push_back(p);
	}

	int szx = myPool.size();
	pos = myPool.pos();

	for (int k = pos; k > 1; --k)
	{
		double* pback = pAlloc[k - 2];
		//pAlloc.pop_back();
		myPool.free(pback);

		pos = myPool.pos();

	}

	//all ok
	//add one
	auto newOne = myPool.alloc();
	(*newOne) = 88;

	myPool.free(newOne);
	myPool.free(pAlloc[MAX_EL - 1]);

	//all deleted
	auto newOnetoo = myPool.alloc();
	(*newOnetoo) = 99;
	//one element 99


	for (int i = 0; i < 3; ++i)
	{
		auto newOnetoo = myPool.alloc();
		(*newOnetoo) = 44 + i;
	}


	//needto test some vakues
}


using namespace DRC::VecD4D;

TEST(TestCaseAlloc, monkyBusinessBuffer) {
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);


	std::vector<double> mix(21,1.0);
	VecXX Vec2(mix);


	auto d = Vec2;
	auto a = d;
	auto b = a;
	auto c = b;


	a *= -1.0;
	auto w = log(-a);
	std::vector<double> cach(w.size());
	for (size_t i = 0; i < w.size(); i++)
	{
		cach[i] = w[i];
	}
	auto aa = -b;

	//operation above should not change
	for (size_t i = 0; i < w.size(); i++)
	{
		double cacI = cach[i];
		double wI = w[i];
		EXPECT_EQ(cacI, wI);
	}

//looks like we need the minus unitary

}