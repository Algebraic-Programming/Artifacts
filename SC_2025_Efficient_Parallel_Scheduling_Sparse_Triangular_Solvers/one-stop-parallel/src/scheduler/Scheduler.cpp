/*
Copyright 2024 Huawei Technologies Co., Ltd.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

@author Toni Boehnlein, Benjamin Lozes, Pal Andras Papp, Raphael S. Steiner   
*/

#include <stdexcept>

#include "scheduler/Scheduler.hpp"

std::pair<RETURN_STATUS, BspSchedule> Scheduler::computeScheduleWithTimeLimit(const BspInstance &instance) {

    std::packaged_task<std::pair<RETURN_STATUS, BspSchedule>(const BspInstance &)> task(
        [this](const BspInstance &instance) -> std::pair<RETURN_STATUS, BspSchedule> {
            return computeSchedule(instance);
        });
    auto future = task.get_future();
    std::thread thr(std::move(task), std::ref(instance));
    if (future.wait_for(std::chrono::seconds(getTimeLimitSeconds())) == std::future_status::timeout) {
        thr.detach(); // we leave the thread still running
        std::cerr << "Timeout reached, execution of computeSchedule() aborted" << std::endl;
        return std::make_pair(TIMEOUT, BspSchedule());
    }
    thr.join();
    try {
        const auto result = future.get();
        return result;
    } catch (const std::exception &e) {
        std::cerr << "Exception caught in computeScheduleWithTimeLimit(): " << e.what() << std::endl;
        return std::make_pair(ERROR, BspSchedule());
    }
}
