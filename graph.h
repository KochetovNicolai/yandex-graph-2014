#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <set>
#include <iostream>

namespace graph {

    template <typename T>
    struct Edge {
        size_t source;
        size_t target;
        T weight;
        Edge() {}
        Edge(size_t source, size_t target, size_t weight);
        Edge backEdge() const;
    };

    template <typename T>
    std::istream &operator>>(std::istream &in, Edge<T> &edge); 
    template <typename T>
    std::ostream &operator<<(std::ostream &out, const Edge<T> &edge); 

    template <typename T>
    struct EdgeList {
        std::vector<Edge<T> > edges;
        // index of the back edge in the target adjacency list
        std::vector<size_t> backEdges;
        Edge<T> &operator[](size_t index);
        const Edge<T> &operator[](size_t index) const;
        size_t size() const;
        void insertEdge(size_t backEdge, const Edge<T> &edge);
        void removeBack();
        size_t findShortestEdge(size_t source, size_t target) const;
    };

    template <typename T>
    class Graph {
        std::vector<EdgeList<T> > adjacencyList;
        // also changes index in the back edge
        void moveEdgeInList(size_t vertex, size_t sourceIndex, size_t targetIndex);
    public:
        Graph() {}
        Graph(size_t size);
        const EdgeList<T> &operator[](size_t index) const;
        size_t size() const;
        void insertEdge(const Edge<T> &edge);
        void removeEdje(size_t vertex, size_t index);
    };

    template <typename T>
    std::istream &operator>>(std::istream &in, Graph<T> &graph);
    template <typename T>
    std::ostream &operator<<(std::ostream &out, const Graph<T> &graph);

    // Dijkstra's algorithm O((V+E)log(V)) for G = <V,E>
    template <typename T, typename Iterator>
    void distanceToSet(const Graph<T> &graph, std::vector<bool> &isAchieved,
                       std::vector<T> &distance, std::vector<size_t> &parent, 
                       Iterator begin, Iterator end);

    // the shortest path between two sets of vertices
    // some path between vertex in source and vertex in target sets
    template <typename T, typename sourceIterator, typename targetIterator>
    T shortestPath(const Graph<T> &graph, std::vector<size_t> &path,
                   sourceIterator sourceBegin, sourceIterator sourceEnd,
                   targetIterator targetBegin, targetIterator targetEnd);

    // Floyd
    template <typename T>
    void pairwiseDistances(const Graph<T> &graph, 
                           std::vector<std::vector<bool> > &isAchieved,
                           std::vector<std::vector<T> > &distance);
}

namespace graph {
    
    template <typename T>
    Edge<T>::Edge(size_t source, size_t target, size_t weight):
        source(source), target(target), weight(weight) {}

    template <typename T>
    Edge<T> Edge<T>::backEdge() const {
        return Edge<T>(target, source, weight);
    }

    template <typename T>
    Edge<T> &EdgeList<T>::operator[](size_t index) { 
        return edges[index];
    }

    template <typename T>
    const Edge<T> &EdgeList<T>::operator[](size_t index) const { 
        return edges[index];
    }

    template <typename T>
    size_t EdgeList<T>::size() const {
        return edges.size();
    }

    template <typename T>
    void EdgeList<T>::insertEdge(size_t backEdge, const Edge<T> &edge) {
        edges.push_back(edge);
        backEdges.push_back(backEdge);
    }

    template <typename T>
    void EdgeList<T>::removeBack() {
        edges.pop_back();
        backEdges.pop_back();
    }

    template <typename T>
    size_t EdgeList<T>::findShortestEdge(size_t source, size_t target) const {
        size_t index = edges.size();
        for (size_t i = 0; i < edges.size(); i++)
            if (edges[i].source == source && edges[i].target == target)
                if (index > i || edges[index].weight > edges[i].weight)
                    index = i;
        return index;
    }

    template <typename T>
    std::istream &operator>>(std::istream &in, Edge<T> &edge) {
        size_t source;
        size_t target;
        size_t weight;
        in >> source >> target >> weight;
        edge = Edge<T>(source, target, weight);
        return in;
    }

    template <typename T>
    std::ostream &operator<<(std::ostream &out, const Edge<T> &edge) {
        out << edge.source << ' ' << edge.target << ' ' << edge.weight;
        return out;
    }

    template <typename T>
    Graph<T>::Graph(size_t size): adjacencyList(size) {}
    
    template <typename T>
    const EdgeList<T> &Graph<T>::operator[](size_t index) const {
        return adjacencyList[index];
    }

    template <typename T>
    size_t Graph<T>::size() const {
        return adjacencyList.size();
    }

    template <typename T>
    void Graph<T>::insertEdge(const Edge<T> &edge) {
        size_t sourceIndex = adjacencyList[edge.source].size();
        size_t targetIndex = adjacencyList[edge.target].size();
        adjacencyList[edge.source].insertEdge(targetIndex, edge);
        adjacencyList[edge.target].insertEdge(sourceIndex, edge.backEdge());
    }

    template <typename T>
    void Graph<T>::moveEdgeInList(size_t vertex, size_t sourceIndex, size_t targetIndex) {
        EdgeList<T> &edgeList = adjacencyList[vertex];
        size_t backEdge = edgeList.backEdges[sourceIndex];
        size_t target = edgeList[sourceIndex].target;
        adjacencyList[target].backEdges[backEdge] = targetIndex;
        edgeList[targetIndex] = edgeList[sourceIndex];
        edgeList.backEdges[targetIndex] = edgeList.backEdges[sourceIndex];
    }

    template <typename T>
    void Graph<T>::removeEdje(size_t vertex, size_t index) {
        size_t target = adjacencyList[vertex][index].target;
        size_t backEdge = adjacencyList[vertex].backEdges[index];
        moveEdgeInList(vertex, adjacencyList[vertex].size() - 1, index);
        moveEdgeInList(target, adjacencyList[target].size() - 1, backEdge);
        adjacencyList[vertex].removeBack();
        adjacencyList[target].removeBack();
    }

    template <typename T>
    std::istream &operator>>(std::istream &in, Graph<T> &graph) {
        size_t size;
        in >> size;
        graph = Graph<T>(size);
        size_t edgesNumber;
        in >> edgesNumber;
        for (size_t i = 0; i < edgesNumber; i++) {
            Edge<T> edge;
            in >> edge;
            graph.insertEdge(edge);
        }
        return in;
    }

    template <typename T>
    std::ostream &operator<<(std::ostream &out, const Graph<T> &graph) {
        out << graph.size();
        for (size_t i = 0; i < graph.size(); i++) {
            for (size_t j = 0; j < graph[i].size(); j++)
                out << graph[i][j];
        }
        return out;
    }

    template <typename T, typename Iterator>
    void distanceToSet(const Graph<T> &graph, std::vector<bool> &isAchieved,
                       std::vector<T> &distance, std::vector<size_t> &parent,
                       Iterator begin, Iterator end) {
        distance.assign(graph.size(), 0);
        parent.assign(graph.size(), 0);
        isAchieved.assign(graph.size(), false);
        std::set< std::pair<size_t, size_t> > queue;
        for (; begin != end; ++begin) {
            size_t vertex = *begin;
            if (!isAchieved[vertex]) {
                parent[vertex] = vertex;
                isAchieved[vertex] = true;
                queue.insert(std::make_pair(0, vertex));
            }
        }
        while (!queue.empty()) {
            size_t vertex = queue.begin()->second;
            queue.erase(queue.begin());
            for (size_t i = 0; i < graph[vertex].size(); i++) {
                Edge<T> edge = graph[vertex][i];
                if (!isAchieved[edge.target] 
                || distance[edge.target] > distance[vertex] + edge.weight) {
                    isAchieved[edge.target] = true;
                    queue.erase(std::make_pair(distance[edge.target], edge.target));
                    distance[edge.target] = distance[vertex] + edge.weight;
                    parent[edge.target] = vertex;
                    queue.insert(std::make_pair(distance[edge.target], edge.target));
                }
            }
        }
    }

    template <typename T, typename sourceIterator, typename targetIterator>
    T shortestPath(const Graph<T> &graph, std::vector<size_t> &path,
                   sourceIterator sourceBegin, sourceIterator sourceEnd,
                   targetIterator targetBegin, targetIterator targetEnd) {
        std::vector<size_t> distance;
        std::vector<size_t> parent;
        std::vector<bool> isAchieved;
        distanceToSet(graph, isAchieved, distance, parent, sourceBegin, sourceEnd);
        bool isFound = false;
        size_t nearestVertex;
        for (; targetBegin != targetEnd; ++targetBegin) {
            size_t vertex = *targetBegin;
            if (!isAchieved[vertex])
                continue;
            if (!isFound || distance[vertex] < distance[nearestVertex])
                nearestVertex = vertex;
            isFound = true;
        }
        path.resize(0);
        if (!isFound)
            return -1;
        while (nearestVertex != parent[nearestVertex]) {
            path.push_back(nearestVertex);
            nearestVertex = parent[nearestVertex];
        }
        path.push_back(nearestVertex);
        for (size_t i = 0; i < path.size() / 2; i++)
            std::swap(path[i], path[path.size() - i - 1]);
        return distance[path.back()];
    }

    template <typename T>
    void pairwiseDistances(const Graph<T> &graph,
                           std::vector<std::vector<bool> > &isAchieved,
                           std::vector<std::vector<T> > &distance) {
        distance.assign(graph.size(), std::vector<T>(graph.size(), T()));
        isAchieved.assign(graph.size(), std::vector<bool>(graph.size(), false));
        for (size_t i = 0; i < graph.size(); i++) {
            for (size_t j = 0; j < graph[i].size(); j++) {
                const Edge<T> &edge = graph[i][j];
                if (!isAchieved[edge.source][edge.target]
                || distance[edge.source][edge.target] > edge.weight) {
                    distance[edge.source][edge.target] = edge.weight;
                    isAchieved[edge.source][edge.target] = true;
                }
            }
        }
        for (size_t i = 0; i < graph.size(); i++) 
            for (size_t j = 0; j < graph.size(); j++) 
                for (size_t k = 0; k < graph.size(); k++) 
                    if (i != j && i != k && isAchieved[j][i] && isAchieved[i][k])
                        if (!isAchieved[j][k] 
                        || distance[j][k] > distance[j][i] + distance[i][k]) {
                              distance[j][k] = distance[j][i] + distance[i][k];
                              isAchieved[j][k] = true;
                        }
    }
}

#endif // GRAPH_H