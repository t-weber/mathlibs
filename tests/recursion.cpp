/**
 * simulating recursion
 * @author Tobias Weber
 * @date December 2021
 * @license see 'LICENSE' file
 *
 * g++ -std=c++20 -Wall -Wextra -Weffc++ -o recursion recursion.cpp
 */

#include <stack>
#include <iostream>
#include <iomanip>


using t_arg = int;

std::stack<t_arg> args;
std::stack<t_arg> rets;

std::stack<t_arg> stack;


// using separate argument and return stacks
void fac_sim()
{
	while(1)
	{
		t_arg arg = args.top(); args.pop();
		t_arg ret = rets.top(); rets.pop();

		if(arg == 0)
			arg = 1;

		// calculation
		rets.push(arg * ret);

		// recursion end condition
		if(arg == 1 || arg == 0)
			break;

		// recursion
		args.push(arg-1);
	}
}


// using a single stack
void fac_sim_2()
{
	while(1)
	{
		t_arg arg = stack.top(); stack.pop();
		t_arg ret = stack.top(); stack.pop();

		if(arg == 0)
			arg = 1;

		// calculation
		stack.push(arg * ret);

		// recursion end condition
		if(arg == 1 || arg == 0)
			break;

		// recursion
		stack.push(arg-1);
	}
}


// actual recursion
t_arg fac_real(t_arg arg)
{
	// end condition
	if(arg == 1 || arg == 0)
		return 1;

	// recursion
	return arg * fac_real(arg-1);
}


// actual recursion
t_arg fibo_real(t_arg arg, t_arg depth=0)
{
	//for(int i=0; i<depth; ++i) std::cout << "    ";
	//std::cout << "arg: " << arg << std::endl;

	if(arg == 0 || arg == 1)
	{
		//for(int i=0; i<depth; ++i) std::cout << "    ";
		//std::cout << "ret: " << 1 << std::endl;
		return 1;
	}

	t_arg ret = fibo_real(arg-1, depth+1) + fibo_real(arg-2, depth+1);

	//for(int i=0; i<depth; ++i) std::cout << "    ";
	//std::cout << "ret: " << ret << std::endl;
	return ret;
}


int main()
{
	for(t_arg i=0; i<10; ++i)
	{
		args.push(i);
		rets.push(1);
		fac_sim();
		int res = rets.top(); rets.pop();

		stack.push(1);
		stack.push(i);
		fac_sim_2();
		t_arg res_2 = stack.top(); stack.pop();

		t_arg res_real = fac_real(i);
		std::cout << "fac(" << i << ") = "
			<< std::setw(10) << std::left << res
			<< std::setw(10) << std::left << res_2
			<< std::setw(10) << std::left << res_real
			<< std::endl;
	}

	std::cout << std::endl;

	for(t_arg i=0; i<6; ++i)
	{
		t_arg res_real = fibo_real(i);
		std::cout << "fibo(" << i << ") = "
			//<< std::setw(10) << std::left << res
			<< std::setw(10) << std::left << res_real
			<< std::endl;
	}

	std::cout << std::endl;

	return 0;
}
