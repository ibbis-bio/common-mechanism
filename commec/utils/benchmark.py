#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Set of tools that allow for a lightweight benchmarking system at the function, and sub-function level.

Usage:
    Create a benchmark_scope("name") instance within any scope you wish to measure.
    The benched time is from instantiation to instance deletion - via either out of scope or manual del.
    Apply an additional stack offset to nested benched scoped, for cleaner visualisation.
    Apply the @benchmark decorator to any function you wish to measure, stack depth handled automatically.

    Note stack depth may collide in visualisation if a nested scope contains @benchmark decorators.
"""
import time
import inspect
import functools

class BenchmarkTool:
    """A simple logger class to handle logging of benchmarks."""
    def __init__(self, filename="benchmark.txt"):
        self.filename = filename
        self.loglines = []
        self.logging_start = time.time()
        self.is_logging = True

    def log(self, message):
        """ Send a log, used by benchmark"""
        if self.is_logging:
            self.loglines.append(message)

    def write_to_file(self):
        """ Write all logs to a file. """
        if self.is_logging:
            self.loglines.reverse()
            with open(self.filename, "w", encoding = 'utf-8"') as file:
                for message in self.loglines:
                    file.write(message + "\n")
        return self.filename

def format_time(seconds):
    """ Simple time formatter, 00:00:00.000 for HH::MM:SS.MS """
    hours, remainder = divmod(seconds, 3600)
    minutes, remainder = divmod(remainder, 60)
    seconds, milliseconds = divmod(remainder, 1)
    milliseconds = int(milliseconds * 1000)  # Convert fraction to milliseconds
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}.{milliseconds:03}"

# This script uses a global single logger instance.
BENCHMARKTOOL = BenchmarkTool()

# Accessors for import to mutate logger.
def benchmark_set_log_file_name(filename : str):
    BENCHMARKTOOL.filename = filename

def benchmark_set_logging(is_logging : bool = True):
    BENCHMARKTOOL.is_logging = is_logging

def benchmark_write_log():
    return BENCHMARKTOOL.write_to_file()

def benchmark_scope(name : str, stack : int = 0):
    """ Factory function for BenchmarkScope instantiation."""
    return BenchmarkScope(name, stack)

class BenchmarkScope:
    """ 
    Benchmarks a given scope from instantiation
    to out of scope or deletion.
    """
    def __init__(self, name : str, stack : int = 0):
        """ 
            Implement a benchmark timer with a given name, 
            and optional additional stack number, 
            for nested scoped benchmarks.
        """
        self.name = name
        self.start_time = time.time() - BENCHMARKTOOL.logging_start
        self.stack_depth : int = int((len(inspect.stack(0))) / 2) + stack

    def __del__(self):
        end_time = time.time() - BENCHMARKTOOL.logging_start
        duration = (end_time - self.start_time)
        log_message = f"{self.stack_depth}\t{self.name}\t{format_time(duration)}\t{format_time(self.start_time)}\t{format_time(end_time)}"
        BENCHMARKTOOL.log(log_message)

# Decorator.
def benchmark(func):
    """
    Decorator to time a function and log its name, stack depth, and timing statistics.
    Writes all information to a log after the program is run. Usefull for gathering statistics
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Measure the start time
        start_time = time.time() - BENCHMARKTOOL.logging_start

        # Get call stack depth
        stack_depth : int = int((len(inspect.stack(0))) / 2)

        # Call the wrapped function
        result = func(*args, **kwargs)

        # Measure the end time and calculate duration
        end_time = time.time() - BENCHMARKTOOL.logging_start
        duration = (end_time - start_time)# * 1000  # Convert to milliseconds

        # Log the function name, depth, and duration
        log_message = f"{stack_depth}\t{func.__name__}\t{format_time(duration)}\t{format_time(start_time)}\t{format_time(end_time)}"
        BENCHMARKTOOL.log(log_message)

        return result
    return wrapper

# Example use case.
if __name__ == "__main__":
    # from utils.benchmark import benchmark, benchmark_set_log_file_name
    benchmark_set_log_file_name("benchmarking.txt")
    @benchmark
    def func1():
        time.sleep(0.1)  # Simulate processing
        func2()
        func2()

    @benchmark
    def func2():
        _bm = benchmark_scope("sleeping")
        time.sleep(1.0)  # Simulate processing
        func3()

    @benchmark
    def func3():
        time.sleep(0.6)  # Simulate processing

    func1()
    benchmark_write_log()
