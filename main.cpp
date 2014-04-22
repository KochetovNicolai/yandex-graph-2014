#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE graphTests

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>

#include "graph.h"
using namespace graph;

size_t checkPath(const Graph<size_t> &graph, std::vector< size_t > &path) {
    if (path.size() == 0)
        return -1;
    for (size_t i = 0; i < path.size(); i++)
        BOOST_CHECK_MESSAGE(path[i] < graph.size(), "nonexistent vertex in path: " << path[i]);
    size_t distance = 0;
    for (size_t i = 1; i < path.size(); i++) {
        size_t index = graph[path[i - 1]].findShortestEdge(path[i - 1], path[i]);
        BOOST_CHECK_MESSAGE(index < graph[path[i - 1]].size(), "nonexistent edge in path: (" 
            << path[i - 1] << ", " << path[i] << ")"); 
        distance += graph[path[i - 1]][index].weight;
    }
    return distance;
}

BOOST_AUTO_TEST_CASE (manualShortestPathTests) {
    std::ifstream in;
    for (size_t i = 0; i < 100; i++) {
        std::string fileName;
        fileName = (std::string("tests/") + char('0' + i / 10)) + char('0' + i % 10);
        in = std::ifstream(fileName + ".in");
        if (!in.is_open()) 
            break;
        Graph<size_t> graph;
        in >> graph;
        size_t sourceSize;
        size_t targetSize;
        in >> sourceSize;
        std::vector<size_t> sourceSet(sourceSize);
        for (size_t i = 0; i < sourceSize; i++)
            in >> sourceSet[i];
        in >> targetSize;
        std::vector<size_t> targetSet(targetSize);
        for (size_t i = 0; i < targetSize; i++)
            in >> targetSet[i];
        std::vector<size_t> path;
        size_t  foundDistance = shortestPath(graph, path,
            sourceSet.begin(), sourceSet.end(), targetSet.begin(), targetSet.end());
        std::ifstream out(fileName + ".out");
        size_t distance;
        out >> distance;
        BOOST_CHECK_EQUAL(distance, foundDistance);
        BOOST_CHECK_EQUAL(distance, checkPath(graph, path));
    }
}


void fillRandomGraph(Graph<size_t> &graph, size_t edgeNumber) {
    for (size_t i = 0; i < edgeNumber; i++) {
        size_t source = (rand() | (rand() << 16)) % graph.size();
        size_t target = (rand() | (rand() << 16)) % graph.size();
        size_t weight = rand();
        graph.insertEdge(Edge<size_t>(source, target, weight));
    }
}

void fillRandomSet(size_t minVal, size_t maxVal, std::vector<size_t> &randomSet) {
    for (size_t i = 0; i < randomSet.size(); i++)
        randomSet[i] =  (rand() | (rand() << 16)) % (maxVal - minVal + 1) + minVal;
}

size_t shortestDistanceViaFloyd(const Graph<size_t> &graph,
                                std::vector<size_t> &source,
                                std::vector<size_t> &target) {
    std::vector<std::vector<bool> > isAchievable;
    std::vector<std::vector<size_t> > distance;
    pairwiseDistances(graph, isAchievable, distance);
    size_t shortest = -1;
    for (size_t i = 0; i < source.size(); i++)
        for (size_t j = 0; j < target.size(); j++)
            if (isAchievable[source[i]][target[j]])
                shortest = min(shortest, distance[source[i]][target[j]]);
    return shortest;
}

BOOST_AUTO_TEST_CASE (randomShortestPathTests) {
    for (size_t i = 2; i < 128; i*= 2) {
        for (size_t j = 0; j < 10; j++) {
            size_t edgeNumber = rand() % (i * i);
            Graph<size_t> graph(i);
            fillRandomGraph(graph, edgeNumber);
            std::vector<size_t> source(rand() % (i / 2 + 1));
            std::vector<size_t> target(rand() % (i / 2 + 1));
            fillRandomSet(0, i / 2 - 1, source);
            fillRandomSet(i / 2, i - 1, target);
            std::vector<size_t> path;
            size_t foundDistance = shortestPath(graph, path,
                source.begin(), source.end(), target.begin(), target.end());
            BOOST_CHECK_EQUAL(foundDistance, checkPath(graph, path));
            size_t distance = shortestDistanceViaFloyd(graph, source, target);
            BOOST_CHECK_EQUAL(foundDistance, distance);
        }
    }
}


BOOST_AUTO_TEST_CASE (largeRandomShortestPathTests) {
    for (size_t i = 1024*512; i <= 1024*1024; i*= 2) {
        for (size_t j = 0; j < 5; j++) {
            size_t edgeNumber = (rand() | (rand() << 16)) % (1024*1024);
            Graph<size_t> graph(i);
            fillRandomGraph(graph, edgeNumber);
            std::vector<size_t> source(rand() % (i / 2 + 1));
            std::vector<size_t> target(rand() % (i / 2 + 1));
            fillRandomSet(0, i / 2 - 1, source);
            fillRandomSet(i / 2, i - 1, target);
            std::vector<size_t> path;
            time_t startTime = clock();
            size_t foundDistance = shortestPath(graph, path,
                source.begin(), source.end(), target.begin(), target.end());
            time_t finishTime = clock();
            double solutionTime = -(startTime - finishTime) / (CLOCKS_PER_SEC + .0);
            BOOST_CHECK_EQUAL(foundDistance, checkPath(graph, path));
            //std::cout << "V = " << i << "  E = " << edgeNumber;
            //std::cout << "  time = " << solutionTime << 's' << std::endl;
        }
    }
}
