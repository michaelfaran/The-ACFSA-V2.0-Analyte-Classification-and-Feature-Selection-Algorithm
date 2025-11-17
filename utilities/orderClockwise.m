function [hrCW, vrCW, idxCW] = orderClockwise(hr, vr)
% orderClockwise  Order 2D points clockwise.
%   [hrCW, vrCW, idxCW] = orderClockwise(hr, vr)
%   hr, vr  : equal-length vectors of x- and y-coordinates
%   hrCW,vrCW : inputs reordered in clockwise order
%   idxCW     : indices such that hrCW = hr(idxCW), vrCW = vr(idxCW)
%
% Notes:
% - Clockwise ordering is obtained by sorting angles around the centroid in
%   descending order. Ties (collinear from centroid) are broken by radius.
% - For collinear sets (all points on a line), “clockwise” is undefined;
%   the function will sort by angle/radius but the geometric orientation
%   has no meaning in that case.

    % sanity
    if numel(hr) ~= numel(vr)
        error('hr and vr must have the same number of elements.');
    end
    hr = hr(:); vr = vr(:);
    N = numel(hr);
    if N < 2
        hrCW = hr.'; vrCW = vr.'; idxCW = (1:N).';
        return
    end

    % centroid
    cx = mean(hr);
    cy = mean(vr);

    % angles and radii w.r.t. centroid
    ang = atan2(vr - cy, hr - cx);           % [-pi, pi]
    rad = hypot(hr - cx, vr - cy);

    % sort: descending angle → clockwise; tie-breaker: farther radius first
    [~, idxCW] = sortrows([ang, rad], [-1, -1]);

    % output (row vectors for convenience)
    hrCW = hr(idxCW).';
    vrCW = vr(idxCW).';
end
