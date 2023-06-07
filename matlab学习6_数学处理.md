# Matlab
## 数据处理
### 线性规划：
目标函数和约束条件均为一次方
**Linprog函数**:
[x, fval] = linprog(f,A,b,Aeq,beq,lb,ub)
* 其中f是目标函数，永列向量表示系数
* A为不等式约束条件的变量系数矩阵，b为常数项矩阵
* Aeq和beq为等式约束
* lb和ub为变量的最小和最大值
* x返回变量的取值，fval返回最值
* linprog只能求最小值，可以先取反

### 非线性规划
**fmincon函数**:
[x, val] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
* fun为脚本定义的目标函数
* x0为变量初始值
* nonlcon为非线性约束

### 多目标规划
可以引入正负偏差变量d+和d-，绝对和目标约束，优先因子等。引入正负偏差后能把大于等于号写为等于号。可以转换为非线性规划问题。
也可以利用fgoalattain函数：
[x, fval, attainfactor, exitflag] = fgoalattain(fun,x0,goal,weight,A,b,Aeq,beq,lb,ub)
* fun为目标函数
* goal为目标，weight为权重

### 最短路径问题：
先建立图和邻接矩阵DG。
sparse函数：DG = sparse([1], [2], W)
* 其中矩阵1和矩阵2对应位置表示对应连接的节点，W表示每条边的权值。
graphshortestpath函数:
[dist, path, pred] = graphshortestpath(DG,起点,终点)
* dist表示最短路径值，path表示路径
```
h = view(biograph(DG,point_name,'ShowWeights','on'))
//绘制图
edges = getedgesbynodeid(h,get(h.Nodes(path),'ID'))
//得到最短路径的边
set(edges,'LineColor',[1 0 0])
```

### 最小生成树
所有节点相互连通
```
//s, t对应的位置为对应节点
G = graph(s,t,weights);
T = minspantree(G)
```

## 路径规划算法：
**1.Dijkstra算法**:
流程：
1. 代价函数g(n)，把代价最低的节点n弹出，标记n为已探索的节点
2. 探索邻居m节点:如果g(m)=infinite，g(m) = g(n) + C~nm~，把m放入队列。如果g(m)>g(n) + C~nm~，令g(m) = g(n) + C~nm~
```
function [ans, next] = matlabdijkstra_3d(adj,begin,destination)
%% adj 所有节点构成的一个可达矩阵，这里假设节点为三维空间中的点
%% begin 源节点序列号
%% destination  目的节点序列号
%% ans  源节点到目的节点的路径长度

% 初始化
U=1:length(adj);  %U存储所有的节点
U(U==begin)=[]; %弹出第一个节点
S=begin;  
dis=adj(begin,:);  %dis表示当前节点到其它节点的距离
curdis=0;  %表示起始节点到当前节点的距离
lastnodes=ones(1,length(dis));  %表示当前节点的前一个节点的数组
%例如lastnodes(5)=4表示5节点前一个是4节点
lastnodes=lastnodes*begin;
%next当前节点最短路径中的下一个节点，如果节点不在最短路径中则为-1
%例如next(3)=2表示节点3的下一个是2
next = -1 * ones(1,length(dis)); % 初始化next数组为-1
while length(S)<length(adj)
    mindis=min(dis(U));  %找到最短的路径
    if mindis==inf
        disp('未完成');
        break
    end
    curdis=curdis+mindis;
    [~,index]=find(dis==mindis,1);  %找到最短路径的节点
    S=[S,index];
    U(U==index)=[];  %弹出该节点
    lastpath=lastnodes(index);  %lastpath记录了找到的节点所经过的上一个节点
    if index==destination
        break;
    end
    
    % 更新最短路径
    for k=1:length(U)
        newdis = curdis + norm(adj(index,:) - adj(U(k),:)); % 3D空间中两点间距离计算
        if(newdis<dis(U(k)))
        %相当于g(m)>g(n)+C~nm~
            dis(U(k))=newdis;
            lastnodes(U(k))=index;
            next(index) = U(k); % 更新next数组
            
        end
    end
end

% 回溯路径
ans=0;
% ~=表示不等于
while destination~=begin
    str=[num2str(next(destination)),'到',num2str(destination)];
    disp(str);
    if next(destination)<0
        disp('找不到合适路径');
        break;
    end
    ans=ans + norm(adj(next(destination),:) - adj(destination,:)); % 计算总路径长度
    destination=next(destination);
end

end
```
这里的g(n)表示两点间距离，利用next数组存放路径并进行回溯。

**2.RRT算法**：
步骤：
1. 随机采样得到新节点rnd
2. 在节点列表node_list中找到与rnd最近的节点nearestnode
3. 计算nearestNode到rnd的方向theta，并生成新节点newNode
4. 检查nearestNode到newNode路径是否与障碍物相交
5. 如果不相交，将newNode加入node_list，检查是否接近目标点，如果是，检查new_node到目标点的路径是否与障碍物相交
6. 如果不相交，返回路径
```
classdef RRT
    properties
        start
        goal
        min_rand
        max_rand
        expand_dis
        goal_sample_rate
        max_iter
        obstacle_list
        node_list
    end
    
    methods
        %初始化RRT函数
        function obj = RRT(obstacleList, randArea, expandDis, goalSampleRate, maxIter)
            obj.start = nan;
            obj.goal = nan;
            %随机采样的最大/最小值
            obj.min_rand = randArea(1);
            obj.max_rand = randArea(2);
            %扩展的步长
            obj.expand_dis = expandDis;
            %目标点被采样的概率
            obj.goal_sample_rate = goalSampleRate;
            obj.max_iter = maxIter;
            obj.obstacle_list = obstacleList;
            obj.node_list = nan;
        end
        
        %rrt路径规划主函数
        function path = rrt_planning(obj,start, goal)
            %将起点和终点节点化，方便计算节点到起点的路径距离
            obj.start = Node(start(1), start(2));
            obj.goal = Node(goal(1), goal(2));
            obj.node_list = [obj.start];
            path = [];
            
            for i = 1:obj.max_iter
                rnd = obj.sample();
                %找到离rnd最近的节点
                n_ind = obj.get_nearset_list_index(obj.node_list, rnd);
                nearestNode = obj.node_list(n_ind);
                
                theta = atan2(rnd(2) - nearestNode.y, rnd(1) - nearestNode.x);
                newNode = obj.get_new_node(theta, n_ind, nearestNode);
                noCollision = obj.check_segment_collision(newNode.x, newNode.y, nearestNode.x, nearestNode.y);
                if noCollision
                    [obj.node_list] = [obj.node_list; newNode];
                    
                    if obj.is_near_goal(newNode)
                        if obj.check_segment_collision(newNode.x, newNode.y, obj.goal.x, obj.goal.y)
                            [rows, ~] = size(obj.node_list);
                            lastIndex = rows;
                            
                            path = obj.get_final_course(lastIndex);
                            pathLen = obj.get_path_len(path);
                            return
                        end
                    end
                end
            end   
        end
        
        function rnd = sample(obj)
        % randi(100)生成1到100之间的随机整数
            if randi(100) > obj.goal_sample_rate
                rnd = [obj.min_rand + (obj.max_rand - obj.min_rand) * rand(), obj.min_rand + (obj.max_rand - obj.min_rand) * rand()];
                % rand生成0到1之间的随机数
            else
                rnd = [obj.goal.x, obj.goal.y];
            end
        end
        
        % 选择新节点的父节点
        function newNode = choose_parent(obj, newNode, nearInds)
            if isempty(nearInds)
                return;
            end
            
            dList = [];
            for i = 1 : length(nearInds)
                % 计算新节点到每个临近节点的距离
                % 如果没有碰撞，则将临近节点的代价加上距离作为该节点的代价
                % 最终选择代价小的节点作为新节点的父节点
                dx = newNode.x - obj.node_list(nearInds(i)).x;
                dy = newNode.y - obj.node_list(nearInds(i)).y;
                d = norm([dx, dy]);
                theta = atan2(dy, dx);
                if obj.check_collision(obj.node_list(i, :), theta, d)
                    dList = [dList; obj.node_list(i, :).cost + d];
                else
                    dList = [dList; inf];
                end
            end
            
            [minCost, minCostIndex] = min(dList);
            minInd = nearInds(minCostIndex);
            
            if minCost == inf
                disp('min cost is inf');
                return
            end
            
            newNode.cost = minCost;
            newNode.parent = minInd;
        end
        

                     
        function Flag = check_segment_collision(obj, x1, y1, x2, y2)
            [rows, ~] = size(obj.obstacle_list);
            for i = 1:rows
                ox = obj.obstacle_list(i, 1); oy = obj.obstacle_list(i, 2); dsize = obj.obstacle_list(i, 3);
                v = [x1, y1];
                w = [x2, y2];
                p = [ox, oy];
                dd = obj.distance_squared_point_to_segment(v, w, p);
                dsize是障碍物大小，相交时v和w的距离小于障碍物大小dsize的平方
                if dd <= dsize^2
                    Flag = false;
                    return
                end
            Flag = true;
            end
        end
        
        function Flag = check_collision(obj, nearNode, theta, d)
            tmpNode = nearNode;
            end_x = tmpNode.x + math.cos(theta) * d;
            end_y = empNode.y + math.sin(theta) * d;
            Flag = obj.check_segment_collision(tmpNode.x, tmpNode.y, end_x, end_y);
        end
        
        % 得到新节点
        function newNode = get_new_node(obj, theta, n_ind, nearestNode)
            newNode = nearestNode;
            newNode.x = newNode.x + obj.expand_dis * cos(theta);
            newNode.y = newNode.y + obj.expand_dis * sin(theta);
            
            newNode.cost = newNode.cost + obj.expand_dis;
            newNode.parent = n_ind;
        end
        
        function Flag = is_near_goal(obj, node)
            d = obj.line_cost(node, obj.goal);
            if d < obj.expand_dis
                Flag = true;
                return
            end
            Flag = false;
        end
        
        function path = get_final_course(obj, lastIndex)
            path = [obj.goal.x, obj.goal.y];
            % 通过parent节点回溯
            while ~isempty(obj.node_list(lastIndex).parent)
                node = obj.node_list(lastIndex);
                path = [path; [node.x, node.y]];
                lastIndex = node.parent;
            end
            path = [path; [obj.start.x, obj.start.y]];
        end
    end
    
    methods(Static)
        function minIndex = get_nearset_list_index(nodes, rnd)
            [rows, ~] = size(nodes);
            dList = [];
            for i = 1:rows
                dList = [dList, (norm([nodes(i).x - rnd(1), nodes(i).y - rnd(2)]))^2];
            end
            [~, minIndex] = min(dList);
        end
        
        function d = line_cost(node1, node2)
            d = norm([node1.x - node2.x, node1.y - node2.y]);
        end
        
        function pathLen = get_path_len(path)
            pathLen = 0;
            [rows, ~] = size(path);
            for i = 2:rows
                node1_x = path(i, 1);
                node1_y = path(i, 2);
                node2_x = path(i-1, 1);
                node2_y = path(i-1, 2);
                pathLen = pathLen + norm([node1_x - node2_x, node1_y - node2_y]);
            end
        end
        
        % 求出p到v，w线段的距离
        % dot表示两向量点积
        function dd = distance_squared_point_to_segment(v, w, p)
            if v == w
                dd = dot(p - v, p - v);
                return
            end
            l2 = dot(w - v, w - v);
            t = max(0, min(1, dot(p - v, w - v) / l2));
            % projection为p在该线段上的投影点
            projection = v + t * (w - v);
            dd = dot(p - projection, p - projection);
        end
    end
end
```

**3.A*算法**：
步骤：
1. 移除f(n) = g(n)+h(n)最小的n节点，标记已发现
2. 对没发现的邻居m节点，采取和dijkstra类似的方法更新g(m)
其中g(n)表示路径的代价，h(n)表示从当前节点到目标节点的估计代价。
注意:
1. h(n)必须要保证小于实际距离，可以采取欧式距离。
2. 可以给h(n)加上权重ω，ω越大表示越趋近于终点方向前进
```

```

**4.人工势场法**：
步骤：
1. 设定起点终点，障碍物，迭代次数，取点半径等
2. 以起点为中心，作半径为r的圆，取n个均布的点
3. 计算n个点前进的代价——终点对其的引力 + 障碍物对其的斥力
4. 取代价最小的点坐标，得到新的起点，重复
5. 发现距离终点很近的点
```
function [ point ] = path_plan(begin,over,obstacle)

iters=1;      % 迭代次数
curr=begin;   % 起点坐标
testR=0.2;    % 测试8点的圆的半径为0.5

while (norm(curr-over)>0.2) &&  (iters<=2000)    % 未到终点&迭代次数不足
   
%     attr=attractive(curr,over);
%     repu=repulsion(curr,obstacle);
    %curoutput=computP(curr,over,obstacle);
    %计算当前点附近半径为0.2的8个点的势能，然后让当前点的势能减去8个点的势能取差值最大的，确定这个
    %方向，就是下一步迭代的点
    
    point(:,iters)=curr;    % point为函数返回值，储存每次遍历得到的新起点 curr
    
    %先求这八个点的坐标
    for i=1:8   % 求 以curr为起点，testR为半径的圆上的八个均匀分布的点
        testPoint(:,i)=[testR*sin((i-1)*pi/4)+curr(1);testR*cos((i-1)*pi/4)+curr(2)];
        testOut(:,i)=computP(testPoint(:,i),over,obstacle);   % 计算上述各个点的所受合力
    end
    [temp num]=min(testOut); % 找出这八个点中，代价最小的点 
    
    %迭代的距离为0.1
    curr=(curr+testPoint(:,num))/2;  % 将上述求得点，迭代到curr上。 （求取的 curr与testPoint 的中点）
    plot(curr(1),curr(2),'og');      % 绘制得到的 新的起点curr
    pause(0.01);            % 程序暂停一会再继续运行 -- 体现出路径搜索的过程
    iters=iters+1;          % 迭代次数+1
end
end


 % 计算周围几个点的势能（代价）
 % 参数：当前起点  终点  障碍物   的坐标
 function [ output ] = computP( curr,over,obstacle )

% 几个用于计算的相关参数 
k_att=1;
repu=0;
k_rep=100;
Q_star=2;     %。障碍物的斥力作用半径

% 计算终点对当前点的引力  
% tips：这个数值越小，说明当前点离终点越近
attr=1/2*k_att*(norm(curr-over))^2;     % 引力计算公式

% 计算所有障碍物对当前点斥力合 
% tips：斥力合越小，说明当前点遇到障碍物的概率越小
for i=1:size(obstacle,2)
    if norm(curr-obstacle(:,i))<=Q_star    % 障碍物到当前点距离在阈值内，考虑其斥力影响
        repu=repu+1/2*k_rep*(1/norm(curr-obstacle(:,i))-1/Q_star)^2;    % 斥力计算公式
        % ps： 当 curr和obstacle 坐标重合时，是否会出现repu计算卡死？ 是否需要对该条件进行设置
    else       % 障碍物到当前点距离在阈值外，忽略斥力影响
        repu=repu+0;
    end
end

output=attr+repu;   % 引力+斥力  这个数值越小，越适合当作下一个起点

 end 
```