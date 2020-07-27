function z=zscore_ignore_nan(r)

zcor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
z=zcor_xnan(r);
