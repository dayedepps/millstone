/**
 * @fileoverview Reference Genome List view.
 */


gd.RefGenomeListView = Backbone.View.extend({
  el: '#gd-page-container',

  initialize: function() {
    this.render();
  },

  render: function() {

    $('#gd-sidenav-link-refgenomes').addClass('active');

    this.datatable = new gd.DataTableComponent({
        el: $('#gd-ref-genome-list-view-datatable-hook'),
        serverTarget: '/_/ref_genomes',
        controlsTemplate: '/_/templates/reference_genome_list_controls',
        requestData: {projectUid: this.model.get('uid')},
    });

    this.uploader = this.$("#uploadDiv").fineUploaderS3({
      debug: true,
      request: {
        endpoint: this.$("#uploadDiv").data("endpoint"),
        accessKey: this.$("#uploadDiv").data("accesskey")
      },
      signature: {
        endpoint: this.$("#uploadDiv").data("signature")
      },
      uploadSuccess: {
        endpoint: this.$("#uploadDiv").data("success")
      },
      objectProperties: {
        key: $.proxy(function(fileId) {
          var filename = this.uploader.fineUploader("getName", fileId);
          return "uploads/" + qq.getUniqueId() + "-" + filename;
        }, this)
      },
      retry: {
        enableAuto: true
      },
      chunking: {
        enabled: true
      },
      deleteFile: {
        endpoint: this.$("#uploadDiv").data("delete"),
        enabled: true,
        forceConfirm: true
      },
      callbacks: {
        onError: function(id, name, reason) {
          alert(reason);
        }
      }
    }).on('complete', $.proxy(function(id, name, response, xhr) {
      var sid = xhr.s3file_id;
      $.post(this.$("#uploadDiv").data("import"), 
            {
              's3file_id': sid,
              'refGenomeLabel': this.$("#uploadDiv").find("#refGenomeLabel").val(),
              'importFileFormat': $("#uploadDiv").find("input[name=importFileFormat]:checked").val(),
            },
            function(data) {
            }
      );}, this));
  },

});
