%STRIPE_FOLD_CHANGE  Calculates the log_2 fold-change in a set of stripe
%values as well as propagating error.
%
%   [FOLDCHANGE FOLDERROR] = stripe_fold_change(STRIPE_MEAN, STRIPE_ERR, CONTROL_MEAN, CONTROL_ERR)
%
%   STRIPE_MEAN and STRIPE_ERR are the mean stripe values and errors
%   respectively.
%   CONTROL_MEAN and CONTROL_ERR are the mean stripe values and errors for
%   the control lines.
%
% Tara Martin, 2014-02-11
% Edited:
%   2014-03-10: changed to include error propagation


function [fold_change fold_error] = stripe_fold_change(stripe_mean, stripe_err, control_mean, control_err)

% first calculate the stripe ratios
[temp temp_err] = errorprop('divide',stripe_mean, stripe_err,... 
                            control_mean, control_err);

% take log_2 of ratios
fold_change = log2(temp);
fold_error = temp_err/(log(2)*temp);

%% OLD VERSION from 2014-02-11
% % get indices of control values
% a = size(values_controls);
% i1 = Icontrols(1);
% i2 = Icontrols(2);
% idx = sub2ind(a, [i1 i1 i2 i2], [1 2 3 4]); %gets the linear indices of a set of rowXcolumn coordinates
% 
% % divide stripe values by control values
% stripes = values(I,:)./values_controls(idx);
% 

