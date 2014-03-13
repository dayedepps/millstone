from django.conf.urls import include
from django.conf.urls import patterns
from django.conf.urls import url
from django.views.generic import RedirectView

import settings

urlpatterns = patterns('',
    url(r'^$', 'main.views.home_view'),

    # Project-specific views
    url(r'^projects$',
            'main.views.project_list_view'),
    url(r'^projects/create$',
            'main.views.project_create_view'),
    url(r'^projects/([\w-]+)$',
            'main.views.project_view'),

    # Tab base views.
    url(r'^projects/([\w-]+)/data$',
            'main.views.project_view'),
    url(r'^projects/([\w-]+)/analyze$',
            'main.views.tab_root_analyze'),
    url(r'^projects/([\w-]+)/analyze/([\w-]+)$',
            'main.views.tab_root_analyze'),
    url(r'^projects/([\w-]+)/analyze/([\w-]+)/([\w-]+)$',
            'main.views.tab_root_analyze'),

    # Reference genomes
    url(r'^projects/([\w-]+)/refgenomes$',
            'main.views.reference_genome_list_view'),
    url(r'^projects/([\w-]+)/refgenomes/([\w-]+)$',
            'main.views.reference_genome_view'),

    # Alignments
    url(r'^projects/([\w-]+)/alignments$',
            'main.views.alignment_list_view'),
    url(r'^projects/([\w-]+)/alignments/create$',
            'main.views.alignment_create_view'),
    url(r'^projects/([\w-]+)/alignments/([\w-]+)$',
            'main.views.alignment_view'),
    url(r'^projects/([\w-]+)/alignments/([\w-]+)/samplealign/([\w-]+)/error$',
            'main.views.sample_alignment_error_view'),

    # Variant sets
    url(r'^projects/([\w-]+)/sets$',
            'main.views.variant_set_list_view'),
    url(r'^projects/([\w-]+)/sets/([\w-]+)$',
            'main.views.variant_set_view'),

    # Samples
    url(r'^projects/([\w-]+)/samples$',
            'main.views.sample_list_view'),

    # Variants
    url(r'^projects/([\w-]+)/refgenomes/([\w-]+)/variants/([\w-]+)$',
            'main.views.single_variant_view'),


    ###########################################################################
    # Templates
    ###########################################################################

    url(r'^templates/sample_list_targets_template.tsv$',
            'main.views.sample_list_targets_template'),

    url(r'^templates/variant_set_upload_template.vcf$',
            'main.views.variant_set_upload_template'),


    ###########################################################################
    # Auth
    ###########################################################################

    # django-registration defaults
    # (further delgates to django.contrib.auth.url)
    (r'^accounts/', include('registration.backends.simple.urls')),

    # The default behavior of registration is redirect to 'users/<username>'.
    # For now let's catch this request here and just redirect to '/'.
    (r'^users/', RedirectView.as_view(url='/')),


    ###########################################################################
    # XHR Actions
    ###########################################################################

    url(r'^_/sets$',
            'main.xhr_handlers.get_variant_set_list'),
    url(r'^_/sets/exportcsv$',
            'main.xhr_handlers.export_variant_set_as_csv'),
    url(r'^_/variants$',
            'main.xhr_handlers.get_variant_list'),
    url(r'^_/variants/modify_set_membership$',
            'main.xhr_handlers.modify_variant_in_set_membership'),
    url(r'^_/variants/refresh_materialized_variant_table$',
            'main.xhr_handlers.refresh_materialized_variant_table'),
    url(r'^_/variants/export_as_csv$',
            'main.xhr_handlers.export_variants_as_csv'),
    url(r'^_/alignmentgroups$',
            'main.xhr_handlers.get_alignment_groups'),
    url(r'^_/samples$',
            'main.xhr_handlers.get_samples'),
    url(r'^_/genes$',
            'main.xhr_handlers.get_gene_list'),
    url(r'^projects/([\w-]+)/delete$',
            'main.xhr_handlers.project_delete'),
    url(r'^_/ref_genomes$',
            'main.xhr_handlers.get_ref_genomes'),


    ###########################################################################
    # Template XHR's
    # TODO: Replace this with client-side templating.
    ###########################################################################

    url(r'^_/templates/variant_filter_controls$',
            'main.template_xhrs.variant_filter_controls'),
    url(r'^_/templates/variant_set_list_controls$',
            'main.template_xhrs.variant_set_list_controls'),
    url(r'^_/templates/alignment_list_controls$',
            'main.template_xhrs.alignment_list_controls'),
    url(r'^_/templates/reference_genome_list_controls$',
            'main.template_xhrs.reference_genome_list_controls'),
    url(r'^_/templates/sample_list_controls$',
            'main.template_xhrs.sample_list_controls'),

)

if settings.S3_ENABLED:
    urlpatterns += patterns('',
        url(r'^_/projects/([\w-]+)/refgenomes/import_s3$',
                'main.xhr_handlers.import_reference_genome_s3',
                name="import_reference_genome_s3"),
        url(r'^_/projects/([\w-]+)/samples/parse_targets_file_s3$',
                'main.xhr_handlers.parse_targets_file_s3',
                name="parse_targets_file_s3"),
        url(r'^_/projects/([\w-]+)/samples/process_sample_files_s3$',
                'main.xhr_handlers.process_sample_files_s3',
                name="process_sample_files_s3"),
        url(r'^s3/signature', 'main.xhr_uploader.handle_s3',
                name="s3_signature"),
        url(r'^s3/delete', 'main.xhr_uploader.handle_s3',
                name='s3_delete'),
        url(r'^s3/success', 'main.xhr_uploader.success',
                name="s3_success")
    )

if settings.RUNNING_ON_EC2:
    urlpatterns += patterns('',
        url(r'^ec2/info$', 'main.views.ec2_info_view',
                name="ec2_info")
    )

if settings.DEBUG:
    from django.conf.urls.static import static
    urlpatterns += static('jbrowse', document_root=settings.JBROWSE_ROOT,
            show_indexes=True)
