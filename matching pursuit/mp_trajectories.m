function traj_mp = mp_trajectories(trajectories)
% Max projection of the trajectories into z dimention
%
% Amin Nejat

    traj_mp = trajectories(:, :, :);
    for t = 1: size(trajectories, 2)
        traj_mp(:, t, 3) = t;
    end

    traj_mp = reshape(traj_mp, numel(traj_mp)/3, []);
    traj_mp(:, 4) = repmat(1:size(trajectories, 1), 1, size(trajectories, 2));
    z_coord = trajectories(:,:,3);
    traj_mp(:, 5) = z_coord(:);
end