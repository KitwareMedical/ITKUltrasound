%define DECL_PYTHON_CurvilinearArraySpecialCoordinatesImage_CLASS(swig_name, template_params)
    %inline %{
    #include "itkContinuousIndexSwigInterface.h"
    %}

    %rename(__SetDirection_orig__) swig_name::SetDirection;
    %extend swig_name {
        itkIndex##template_params TransformPhysicalPointToIndex( itkPointD##template_params & point ) {
            itkIndex##template_params idx;
            self->TransformPhysicalPointToIndex<double>( point, idx );
            return idx;
         }

        itkContinuousIndexD##template_params TransformPhysicalPointToContinuousIndex( itkPointD##template_params & point ) {
            itkContinuousIndexD##template_params idx;
            self->TransformPhysicalPointToContinuousIndex<double>( point, idx );
            return idx;
        }

        itkPointD##template_params TransformContinuousIndexToPhysicalPoint( itkContinuousIndexD##template_params & idx ) {
            itkPointD##template_params point;
            self->TransformContinuousIndexToPhysicalPoint<double>( idx, point );
            return point;
        }

        itkPointD##template_params TransformIndexToPhysicalPoint( itkIndex##template_params & idx ) {
            itkPointD##template_params point;
            self->TransformIndexToPhysicalPoint<double>( idx, point );
            return point;
        }
   }

%enddef
