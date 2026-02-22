// Copyright Gradientspace Corp. All Rights Reserved.
#pragma once

#include "Core/gs_parallel_api.h"

namespace GS
{

class DummyTaskWrapper : public IExternalTaskWrapper
{
public:
	DummyTaskWrapper() {}
};

class DummyParallelAPI : public Parallel::gs_parallel_api
{
public:
	void parallel_for_jobcount(
		uint32_t NumJobs,
		FunctionRef<void(uint32_t JobIndex)> JobFunction,
		ParallelForFlags Flags) override
	{
		for (uint32_t i = 0; i < NumJobs; i++)
		{
			JobFunction(i);
		}
	}

	TaskContainer launch_task(
		const char* Identifier,
		std::function<void()> task,
		TaskFlags Flags) override
	{
		task();
		TaskContainer Container;
		Container.ExternalTask = std::make_shared<DummyTaskWrapper>();
		return Container;
	}

	void wait_for_task(TaskContainer& task) override
	{
		// task already completed in launch_task
	}
};

}
