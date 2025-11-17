function[average_distance]= calc_avg_dist(points)
% % Define the coordinates of the points
% points = [
%     1, 1;
%     -1, 1;
% 
% ];

% Number of points
num_points = size(points, 1);

% Initialize the sum of distances
sum_distances = 0;

% Calculate pairwise distances
for i = 1:num_points
    for j = i+1:num_points
        % Calculate the distance between point i and point j
        distance = sqrt((points(i, 1) - points(j, 1))^2 + (points(i, 2) - points(j, 2))^2);
        % Add to the sum of distances
        sum_distances = sum_distances + distance;
    end
end

% Calculate the number of unique pairs
num_pairs = nchoosek(num_points, 2);

% Calculate the average distance
average_distance = sum_distances / num_pairs;

% Display the result
% fprintf('The average distance between all points is: %.4f\n', average_distance);
