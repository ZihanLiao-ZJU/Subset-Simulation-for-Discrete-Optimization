function view = ISDIntegrationView(intBay)
% Density evaluator in the measure used for resampling and integration.
% Identity by default; projected ISDs supply an IntegrationView value object.
space=ISDSpace(intBay);
view=space.IntegrationView(intBay);
end
